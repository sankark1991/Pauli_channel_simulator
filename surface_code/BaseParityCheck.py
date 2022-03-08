import numpy as np
from copy import copy

def totuple(a):
    try: return tuple(totuple(i) for i in a)
    except TypeError: return a

def toarray(a):
    try: return np.stack(toarray(i) for i in a)
    except TypeError: return a

class BaseParityCheck:
    """Rotated surface code with distance 5."""
    def __init__(self):
        """
        Generates the parity check circuits.
        Gets and stores the full list of faults by type, round number, qubits, and Pauli.
        Pre-computes and stores the error and syndrome difference for each one.
        For each edge (v1, v2), stores the list of faults (l, P) which trigger that edge."""
        self.get_parity_check_circuit()
        print("Generated parity check circuit")
        self.get_parity_check_circuit_array_form()
        print("Converted parity check circuit to array form")

        self.get_all_fault_locations()
        print("Generated list of fault locations")
        print(f"Number of R faults: {len(self.R_faults)}")
        print(f"Number of M faults: {len(self.M_faults)}")
        print(f"Number of MR faults: {len(self.MR_faults)}")
        print(f"Number of CNOT_1 faults: {len(self.CNOT_1_faults)}")
        print(f"Number of CNOT_2 faults: {len(self.CNOT_2_faults)}")
        print(f"Number of IR faults: {len(self.IR_faults)}")
        print(f"Number of IM faults: {len(self.IM_faults)}")
        print(f"Number of IMR faults: {len(self.IMR_faults)}")
        print(f"Number of IC faults: {len(self.IC_faults)}")


        self.R_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.R_faults}
        self.M_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.M_faults}
        self.MR_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.MR_faults}
        self.CNOT_1_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.CNOT_1_faults}
        self.CNOT_2_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.CNOT_2_faults}
        self.IR_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.IR_faults}
        self.IM_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.IM_faults}
        self.IMR_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.IMR_faults}
        self.IC_fault_effects = {fault: self.error_and_syn_diff(*fault) for fault in self.IC_faults}
        print("Calculated effect of all possible faults")

        # TODO edge counts

    def get_parity_check_circuit(self):
        """Encodes the parity check circuit, and also stores the list of error locations/types (l, P).
        The CNOTs can be scheduled so that:
            - no syndrome qubit is ever idle
            - data qubits in the bulk are never idle
            - data qubits along the bottom and left edge are idle in round 1
            - data qubits along the top and right edge are idle in round 4
            - data qubit in the bottom left corner is isle in rounds 1 and 2
            - data qubit in the bottom right and top left corners are idle in rounds 1 and 4
            - data qubit in the top right corner is idle in rounds 3 and 4"""
        self.data = [np.array([2 * i + 1, 2 * j + 1]) for j in range(5) for i in range(5)]
        self.X_syns = [np.array([2 + 4 * i - 2 * k, 2 + 4 * j + 2 * k])
                       for k in range(2) for j in range(2) for i in range(3)]
        self.Z_syns = [np.array([2 + 4 * j + 2 * k, 4 * i + 2 * k])
                       for k in range(2) for j in range(2) for i in range(3)]
        # print(self.data)
        # print(self.X_syns)
        # print(self.Z_syns)
        # The part of the parity check circuit where noise can occur
        round_0 = {'RX': self.X_syns, 'RZ': self.Z_syns, 'IR': self.data}

        round_1 = {
            'CNOT': [(v, v + np.array([1, 1])) for v in self.X_syns if 1 < v[0] < 9] + [(v, v + np.array([1, 1])) for v
                                                                                        in
                                                                                        self.Z_syns if 1 < v[1] < 9],
            'IC': [q for q in self.data if q[0] == 1 or q[1] == 1]}

        round_2 = {
            'CNOT': [(v, v + np.array([-1, 1])) for v in self.X_syns if 1 < v[0]] + [(v, v + np.array([1, 1])) for v in
                                                                                     self.X_syns if v[0] < 1] + [
                        (v, v + np.array([1, -1])) for v in
                        self.Z_syns if 1 < v[1]] + [(v, v + np.array([1, 1])) for v in
                                                    self.Z_syns if v[1] < 1],
            'IC': [np.array([1, 1])]}
        round_3 = {
            'CNOT': [(v, v + np.array([1, -1])) for v in self.X_syns if v[0] < 9] + [(v, v + np.array([-1, -1])) for v
                                                                                     in self.X_syns if 9 < v[0]] + [
                        (v, v + np.array([-1, 1])) for v in
                        self.Z_syns if v[1] < 9] + [(v, v + np.array([-1, -1])) for v in
                                                    self.Z_syns if 9 < v[1]],
            'IC': [np.array([9, 9])]}

        round_4 = {
            'CNOT': [(v, v + np.array([-1, -1])) for v in self.X_syns if 1 < v[0] < 9] + [(v, v + np.array([-1, -1]))
                                                                                          for v in
                                                                                          self.Z_syns if 1 < v[1] < 9],
            'IC': [q for q in self.data if q[0] == 9 or q[1] == 9]}

        round_5 = {'MRX': self.X_syns, 'MRZ': self.Z_syns, 'IMR': self.data}

        # The part of the parity check circuit which is perfect
        round_6 = {
            'CNOT': [(v, v + np.array([1, 1])) for v in self.X_syns if 1 < v[0] < 9] + [(v, v + np.array([1, 1])) for v
                                                                                        in
                                                                                        self.Z_syns if 1 < v[1] < 9]}

        round_7 = {
            'CNOT': [(v, v + np.array([-1, 1])) for v in self.X_syns if 1 < v[0]] + [(v, v + np.array([1, 1])) for v in
                                                                                     self.X_syns if v[0] < 1] + [
                        (v, v + np.array([1, -1])) for v in
                        self.Z_syns if 1 < v[1]] + [(v, v + np.array([1, 1])) for v in
                                                    self.Z_syns if v[1] < 1]}

        round_8 = {
            'CNOT': [(v, v + np.array([1, -1])) for v in self.X_syns if v[0] < 9] + [(v, v + np.array([-1, -1])) for v
                                                                                     in self.X_syns if 9 < v[0]] + [
                        (v, v + np.array([-1, 1])) for v in
                        self.Z_syns if v[1] < 9] + [(v, v + np.array([-1, -1])) for v in
                                                    self.Z_syns if 9 < v[1]]}

        round_9 = {
            'CNOT': [(v, v + np.array([-1, -1])) for v in self.X_syns if 1 < v[0] < 9] + [(v, v + np.array([-1, -1]))
                                                                                          for v in
                                                                                          self.Z_syns if 1 < v[1] < 9]}

        round_10 = {'MX': self.X_syns, 'MZ': self.Z_syns}

        # Store the circuit and fault locations
        self.circuit = [round_0, round_1, round_2, round_3, round_4, round_5, round_6, round_7, round_8, round_9,
                        round_10]
        self.fault_locations = []
        for round_number in range(6):
            round = self.circuit[round_number]
            for fault_type in round.keys():
                target_qubits_list = round[fault_type]
                for qbs in target_qubits_list:
                    self.fault_locations.append((round_number, fault_type, qbs))
        return None

    def get_parity_check_circuit_array_form(self):
        """Store parity check circuit as an array, with the following conventions:
            - RX/RZ = prepare X/Z
            - MX/MZ = measure X/Z
            - MRX/MRZ = measure and reset X/Z
            - 1000 + n = control of CNOT with target n
            - 2000 + n = target of CNOT with control n
            - IR/IM/IMR/IC = idle
            - 0 = protected location (no noise)"""
        all_qubits = self.data + self.X_syns + self.Z_syns
        # print(all_qubits)
        self.num_qubits = len(all_qubits)
        # Create a lookup array, shifted by 1 to make later operations faster
        min_x, max_x = min(q[0] for q in all_qubits), max(q[0] for q in all_qubits)
        min_y, max_y = min(q[1] for q in all_qubits), max(q[1] for q in all_qubits)
        # print(min_x, max_x)
        # print(min_y, max_y)
        self.qubit_ind_array = np.array([[None for _ in range(max_y - min_y + 1)]
                                         for _ in range(max_x - min_x + 1)])
        for q in all_qubits:
            # print(q)
            q_ind = [i for i in range(self.num_qubits) if np.allclose(q, all_qubits[i])][0]
            self.qubit_ind_array[q[0], q[1]] = q_ind
        # print(self.qubit_ind_array)
        self.data_inds = [self.qubit_ind_array[totuple(q)] for q in self.data]
        self.X_syn_inds = [self.qubit_ind_array[totuple(q)] for q in self.X_syns]
        self.Z_syn_inds = [self.qubit_ind_array[totuple(q)] for q in self.Z_syns]

        self.circuit_array = np.array([[None for _ in range(self.num_qubits)] for _ in range(11)])

        for round_number in range(11):
            round = self.circuit[round_number]
            for op_type in round.keys():
                # print(op_type)
                if any(op_type == ft for ft in ('MX', 'MZ', 'RX', 'RZ', 'MRX', 'MRZ', 'IR', 'IM', 'IMR', 'IC')):
                    for v in round[op_type]:
                        # print(v)
                        # print(self.qubit_ind_array[v[0], v[1]])
                        self.circuit_array[round_number, self.qubit_ind_array[totuple(v)]] = op_type
                elif op_type == 'CNOT':
                    for g in round['CNOT']:
                        # print(round_number, g, tuple(self.qubit_ind_array[totuple(x)] for x in g))
                        self.circuit_array[round_number, self.qubit_ind_array[totuple(g[0])]] = 1000 + \
                                                                                                self.qubit_ind_array[
                                                                                                    totuple(g[1])]
                        self.circuit_array[round_number, self.qubit_ind_array[totuple(g[1])]] = 2000 + \
                                                                                                self.qubit_ind_array[
                                                                                                    totuple(g[0])]
        return None

    def get_all_fault_locations(self):
        """Gets all fault locations"""
        # R: Pauli after a preparation operation
        # M, MR: Pauli before a measurement operations
        # CNOT: 2-Pauli after a CNOT gate. Two possible dists to draw from
        # IR, IM, IC: Pauli after an idling operation
        self.R_faults = []
        self.M_faults = []
        self.MR_faults = []
        self.CNOT_1_faults = []
        self.CNOT_2_faults = []
        self.IR_faults = []
        self.IM_faults = []
        self.IMR_faults = []
        self.IC_faults = []
        for round_num in range(6):
            round = self.circuit[round_num]
            if 'RX' in round.keys():
                self.R_faults.extend([(round_num, (totuple(q),), (P,)) for q in round['RX'] for P in ['Z', 'Y']])
            if 'RZ' in round.keys():
                self.R_faults.extend([(round_num, (totuple(q),), (P,)) for q in round['RZ'] for P in ['X', 'Y']])
            if 'MX' in round.keys():
                self.M_faults.extend([(round_num - 1, (totuple(q),), (P,)) for q in round['MX'] for P in ['Z', 'Y']])
            if 'MZ' in round.keys():
                self.M_faults.extend([(round_num - 1, (totuple(q),), (P,)) for q in round['MZ'] for P in ['X', 'Y']])
            if 'MRX' in round.keys():
                self.MR_faults.extend([(round_num - 1, (totuple(q),), (P,)) for q in round['MRX'] for P in ['Z', 'Y']])
            if 'MRZ' in round.keys():
                self.M_faults.extend([(round_num - 1, (totuple(q),), (P,)) for q in round['MRZ'] for P in ['X', 'Y']])
            if 'IR' in round.keys():
                self.IR_faults.extend([(round_num, (totuple(q),), (P,)) for q in round['IR'] for P in ['Z', 'X', 'Y']])
            if 'IM' in round.keys():
                self.IM_faults.extend([(round_num, (totuple(q),), (P,)) for q in round['IM'] for P in ['Z', 'X', 'Y']])
            if 'IMR' in round.keys():
                self.IMR_faults.extend(
                    [(round_num, (totuple(q),), (P,)) for q in round['IMR'] for P in ['Z', 'X', 'Y']])
            if 'IC' in round.keys():
                self.IC_faults.extend([(round_num, (totuple(q),), (P,)) for q in round['IC'] for P in ['Z', 'X', 'Y']])
            if 'CNOT' in round.keys():
                self.CNOT_1_faults.extend([(round_num, ((totuple(Q[i])),), (P,)) for Q in round['CNOT']
                                           for P in ['Z', 'X', 'Y'] for i in range(2)])
                self.CNOT_2_faults.extend([(round_num, totuple(Q), (P_1, P_2)) for Q in round['CNOT']
                                           for P_1 in ['I', 'Z', 'X', 'Y']
                                           for P_2 in ['I', 'Z', 'X', 'Y'] if not (P_1 == 'I' and P_2 == 'I')])

    def error_and_syn_diff(self, r, Q, P):
        """Given an error (r, Q, P), where r is the round number, Q is the qubits afflicted (in 2D coordinate form),
        P is the Pauli error, returns the net effect on the data qubits and the syndrome diff history.
        Computed by propagating through the base parity check circuit in array form."""
        # Generate the error state after round r.
        pauli_binary = {'I': np.array([0, 0]), 'X': np.array([0, 1]), 'Z': np.array([1, 0]), 'Y': np.array([1, 1])}
        error = np.zeros((self.num_qubits, 2))
        for i in range(len(Q)):
            error[self.qubit_ind_array[totuple(Q[i])]] = pauli_binary[P[i]]

        syn = np.zeros((2, self.num_qubits))
        # Propagate forward until they get measured out
        after_round = r
        while any(not np.allclose(error[i], np.zeros(2)) for i in self.X_syn_inds + self.Z_syn_inds):
            # While there are still errors on any syndrome indices, propagate forward
            after_round += 1  # We are going to propagate past the next round
            qubits_to_handle = [i for i in range(self.num_qubits)]
            while len(qubits_to_handle) > 0:
                i = qubits_to_handle[0]
                # Figure out what gate qubit i is involved in.
                if any(self.circuit_array[after_round, i] == op_type for op_type in ('MX', 'MRX')):
                    # Measure X
                    if r < 5:
                        syn[0, i] = error[i, 0]  # If Z or Y, flips measurement outcome in first syndrome round
                    elif r >= 5:
                        syn[1, i] = error[i, 0]  # If Z or Y, flips measurement outcome in second syndrome round
                    error[i] = np.array([0, 0])
                elif any(self.circuit_array[after_round, i] == op_type for op_type in ('MZ', 'MRZ')):
                    # Measure Z
                    if r < 5:
                        syn[0, i] = error[i, 1]  # If X or Y, flips measurement outcome in first syndrome round
                    elif r >= 5:
                        syn[1, i] = error[i, 1]  # If X or Y, flips measurement outcome in second syndrome round
                    error[i] = np.array([0, 0])
                elif any(self.circuit_array[after_round, i] == op_type for op_type in
                         ('RX', 'RZ', 'IR', 'IM', 'IMR', 'IC', None)):
                    pass
                elif 1000 <= self.circuit_array[after_round, i] < 2000:
                    # CNOT gate from i to j
                    j = self.circuit_array[after_round, i] - 1000
                    error[i] = error[i] + error[j, 0] * np.array([1, 0])
                    error[j] = error[j] + error[i, 1] * np.array([0, 1])
                    qubits_to_handle.remove(j)
                elif self.circuit_array[after_round, i] >= 2000:
                    # CNOT gate from j to i
                    j = self.circuit_array[after_round, i] - 2000
                    error[j] = error[j] + error[i, 0] * np.array([1, 0])
                    error[i] = error[i] + error[j, 1] * np.array([0, 1])
                    qubits_to_handle.remove(j)
                qubits_to_handle.remove(i)
        syn_diff = np.stack((syn[0], np.mod(syn[1] + syn[0], 2)))

        return error, syn_diff

    def edge_error_list(self, e1, e2):
        """Given vertices v1 = (x1, y1, t1) and v2 = (x2, y2, t2) with t1, t2 in {0, 1}, returns the set of
        faults (r, Q, P) which trigger this edge."""
        pass
