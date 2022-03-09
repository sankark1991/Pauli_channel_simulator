"""
Here's the full algorithm.

1. First, define RQP format for faults. RQP format means a 5-tuple
    (round number the fault is after, x-coord of qubit, y-coord of qubit, Z part of Pauli, X part of Pauli)
For d=5, the round number is in [0, 5], the coords are in [0, 10], and the Pauli parts in [0, 1]

2. For the d=5 surface code, pre-compute the full mapping from RQP format to the resulting errors and sdiffs.
The error and sdiff are arrays of shape (2, 11, 11). Store this all in a dictionary, "results_dict"

3. Now, the sampling algorithm:
- Generate the d x d surface code circuit
- Sample a list of faults, based on some noise model with parameters (p1, p2, ...).
- Convert the faults to RQP format. CNOT faults can be broken into two RQP faults.
- Shrink down to 5x5 surface code, get effect, grow back up, add to total

"""
import numpy as np
import random
import pickle
from .BaseParityCheck import BaseParityCheck, totuple, toarray

def shrink_space(v, d):
    """Given code distance d, shrinks coordinates to those appropriate for code distance 5"""
    try: return tuple(shrink_space(n, d) for n in v)
    except TypeError:
        if v <= 3:
            return v
        elif v >= 2 * d - 3:
            return v - 2 * d + 10
        else:
            x = (v - 3) // 4
            return v - 4 * x

def shrink_time(t, N_rounds):
    """Given N_rounds, shrinks time coordinate to those appropriate for a 1 noisy + 1 perfect round circuit"""
    if t == 0:
        return 0
    elif t <= 5 * N_rounds:
        x = (5 - t) // 5
        return t + 5 * x
    else:
        return t + 5 - 5 * N_rounds

class ParityCheck:
    def __init__(self, code_dist, N_rounds, small_code):
        self.code_dist = code_dist
        self.N_rounds = N_rounds
        self.small_code = small_code
        self.construct_circuit()

    def construct_circuit(self, op_types=('RX', 'RZ', 'IR', 'MRX', 'MRZ', 'IMR', 'MX', 'MZ', 'IM', 'CNOT', 'IC')):
        data = [np.array([2 * i + 1, 2 * j + 1]) for j in range(self.code_dist) for i in range(self.code_dist)]
        X_syns = [np.array([2 + 4 * i - 2 * k, 2 + 4 * j + 2 * k])
                  for k in range(2) for j in range((self.code_dist - 1) // 2) for i in range((self.code_dist + 1) // 2)]
        Z_syns = [np.array([2 + 4 * j + 2 * k, 4 * i + 2 * k])
                  for k in range(2) for j in range((self.code_dist - 1) // 2) for i in range((self.code_dist + 1) // 2)]
        round_0 = {'RX': X_syns, 'RZ': Z_syns, 'IR': data}
        round_1 = {
            'CNOT': [(v, v + np.array([1, 1])) for v in X_syns if 1 < v[0] < 2 * self.code_dist - 1] + [(v, v + np.array([1, 1])) for v
                                                                                                        in
                                                                                                        Z_syns if 1 < v[1] < 2 * self.code_dist - 1],
            'IC': [q for q in data if q[0] == 1 or q[1] == 1]}

        round_2 = {
            'CNOT': [(v, v + np.array([-1, 1])) for v in X_syns if 1 < v[0]] + [(v, v + np.array([1, 1])) for v in
                                                                                X_syns if v[0] < 1] + [
                        (v, v + np.array([1, -1])) for v in
                        Z_syns if 1 < v[1]] + [(v, v + np.array([1, 1])) for v in
                                               Z_syns if v[1] < 1],
            'IC': [np.array([1, 1])]}
        round_3 = {
            'CNOT': [(v, v + np.array([1, -1])) for v in X_syns if v[0] < 2 * self.code_dist - 1] + [(v, v + np.array([-1, -1])) for v
                                                                                                     in X_syns if 2 * self.code_dist - 1 < v[0]] + [
                        (v, v + np.array([-1, 1])) for v in
                        Z_syns if v[1] < 2 * self.code_dist - 1] + [(v, v + np.array([-1, -1])) for v in
                                                                    Z_syns if 2 * self.code_dist - 1 < v[1]],
            'IC': [np.array([9, 9])]}

        round_4 = {
            'CNOT': [(v, v + np.array([-1, -1])) for v in X_syns if 1 < v[0] < 2 * self.code_dist - 1] + [(v, v + np.array([-1, -1]))
                                                                                                          for v in
                                                                                                          Z_syns if 1 < v[1] < 2 * self.code_dist - 1],
            'IC': [q for q in data if q[0] == 2 * self.code_dist - 1 or q[1] == 2 * self.code_dist - 1]}

        round_5 = {'MRX': X_syns, 'MRZ': Z_syns, 'IMR': data}

        self.circuit = [round_0, round_1, round_2, round_3, round_4, round_5]
        for round in self.circuit:
            for op_type in op_types:
                if op_type not in round.keys():
                    round[op_type] = []

    def add_translate(self, r, q0, q1, pz, px):
        """Take an RQP for the big code, go down to an RQP for the small code, retrieve results, translate back up and add."""
        # Shrink step
        r_shrink = shrink_time(r, self.N_rounds)
        q0_shrink, q1_shrink = shrink_space((q0, q1), self.code_dist)

        # Look up results.
        error = self.small_code.fault_error[(r_shrink, q0_shrink, q1_shrink, pz, px)]
        sdiff = self.small_code.fault_sdiff[(r_shrink, q0_shrink, q1_shrink, pz, px)]
        # error, sdiff = results_dict[(r_shrink, *Q_shrink, *P)]

        # Translate and add
        r_tr = (r - r_shrink) // 5
        Q_xtr = q0 - q0_shrink
        Q_ytr = q1 - q1_shrink
        self.full_error[Q_xtr: Q_xtr + 11, Q_ytr: Q_ytr + 11] += error
        self.full_sdiff[r_tr: r_tr + 2, Q_xtr: Q_xtr + 11, Q_ytr: Q_ytr + 11] += sdiff

    def sample(self, p):
        error_rates = {'RX': p, 'RZ': p, 'MRX': p, 'MRZ': p, 'MX': p, 'MZ': p, 'CNOT_1': p, 'CNOT_2': p,
                       'IR': p, 'IM': p, 'IMR': p, 'IC': p}
        self.full_error = np.zeros((2 * self.code_dist + 1, 2 * self.code_dist + 1, 2))
        self.full_sdiff = np.zeros((self.N_rounds + 1, 2 * self.code_dist + 1, 2 * self.code_dist + 1))
        RQPs = []
        for round_num in range(6):
            round = self.circuit[round_num]
            for t in range(self.N_rounds):
                r = round_num + 5 * t
                # Preparation errors and associated idling
                for Q in round['RX']:
                    z = random.random()
                    if z < error_rates['RX'] / 4:
                        RQPs.append((r, Q[0], Q[1], 1, 0))
                    elif z < error_rates['RX'] / 2:
                        RQPs.append((r, Q[0], Q[1], 1, 1))
                for Q in round['RZ']:
                    z = random.random()
                    if z < error_rates['RZ'] / 4:
                        RQPs.append((r, Q[0], Q[1], 0, 1))
                    elif z < error_rates['RZ'] / 2:
                        RQPs.append((r, Q[0], Q[1], 1, 1))
                for Q in round['IR']:
                    z = random.random()
                    if z < error_rates['IR'] / 4:
                        RQPs.append((r, Q[0], Q[1], 1, 0))
                    elif z < error_rates['IR'] / 2:
                        RQPs.append((r, Q[0], Q[1], 0, 1))
                    elif z < 3 * error_rates['IR'] / 4:
                        RQPs.append((r, Q[0], Q[1], 1, 1))
                # Measurement-reset errors and associated idling
                for Q in round['MRX']:
                    z = random.random()
                    if z < error_rates['MRX'] / 4:
                        RQPs.append((r-1, Q[0], Q[1], 1, 0))
                    elif z < error_rates['MRX'] / 2:
                        RQPs.append((r-1, Q[0], Q[1], 1, 1))
                for Q in round['MRZ']:
                    z = random.random()
                    if z < error_rates['MRZ'] / 4:
                        RQPs.append((r-1, Q[0], Q[1], 0, 1))
                    elif z < error_rates['MRZ'] / 2:
                        RQPs.append((r-1, Q[0], Q[1], 1, 1))
                for Q in round['IMR']:
                    z = random.random()
                    if z < error_rates['IMR'] / 4:
                        RQPs.append((r-1, Q[0], Q[1], 1, 0))
                    elif z < error_rates['IMR'] / 2:
                        RQPs.append((r-1, Q[0], Q[1], 0, 1))
                    elif z < 3 * error_rates['IMR'] / 4:
                        RQPs.append((r-1, Q[0], Q[1], 1, 1))
                # Measurement errors and associated idling
                for Q in round['MX']:
                    z = random.random()
                    if z < error_rates['MX'] / 4:
                        RQPs.append((r-1, Q[0], Q[1], 1, 0))
                    elif z < error_rates['MX'] / 2:
                        RQPs.append((r-1, Q[0], Q[1], 1, 1))
                for Q in round['MZ']:
                    z = random.random()
                    if z < error_rates['MZ'] / 4:
                        RQPs.append((r-1, Q[0], Q[1], 0, 1))
                    elif z < error_rates['MZ'] / 2:
                        RQPs.append((r-1, Q[0], Q[1], 1, 1))
                for Q in round['IM']:
                    z = random.random()
                    if z < error_rates['IM'] / 4:
                        RQPs.append((r-1, Q[0], Q[1], 1, 0))
                    elif z < error_rates['IM'] / 2:
                        RQPs.append((r-1, Q[0], Q[1], 0, 1))
                    elif z < 3 * error_rates['IM'] / 4:
                        RQPs.append((r-1, Q[0], Q[1], 1, 1))
                # CNOT errors and associated idling
                for Qs in round['CNOT']:
                    for Q in Qs:
                        z = random.random()
                        if z < error_rates['CNOT_1'] / 4:
                            RQPs.append((r, Q[0], Q[1], 1, 0))
                        elif z < error_rates['CNOT_1'] / 2:
                            RQPs.append((r, Q[0], Q[1], 0, 1))
                        elif z < 3 * error_rates['CNOT_1'] / 4:
                            RQPs.append((r, Q[0], Q[1], 1, 1))
                for Q in round['IC']:
                    z = random.random()
                    if z < error_rates['IC'] / 4:
                        RQPs.append((r, Q[0], Q[1], 1, 0))
                    elif z < error_rates['IC'] / 2:
                        RQPs.append((r, Q[0], Q[1], 0, 1))
                    elif z < 3 * error_rates['IC'] / 4:
                        RQPs.append((r, Q[0], Q[1], 1, 1))
        for rqp in RQPs:
            self.add_translate(*rqp)
        self.full_error = np.mod(self.full_error, 2)
        self.full_sdiff = np.mod(self.full_sdiff, 2)

        return self.full_error, self.full_sdiff