import numpy as np

class BaseParityCheck:
    """Rotated surface code with distance 5."""
    def __init__(self):
        """Gets and stores the full list of error locations/types (l, P).
        Pre-computes and stores syn_diff(l, P) and error(l, P) for each one.
        Then for each edge (v1, v2), stores the list of faults (l, P) which trigger that edge."""
        pass

    def get_circuit(self):
        """Encodes the parity check circuit, and also stores the list of error locations/types (l, P).
        The CNOTs can be scheduled so that:
            - no syndrome qubit is ever idle
            - data qubits in the bulk are never idle
            - data qubits along the bottom and left edge are idle in round 1
            - data qubits along the top and right edge are idle in round 4
            - data qubit in the bottom left corner is isle in rounds 1 and 2
            - data qubit in the bottom right and top left corners are idle in rounds 1 and 4
            - data qubit in the top right corner is idle in rounds 3 and 4"""
        self.data = [(2 * i, 2 * j) for j in range(5) for i in range(5)]
        self.X_syns = [np.array([1 + 4 * i - 2 * k, 1 + 4 * j + 2 * k])
                       for k in range(2) for j in range(2) for i in range(3)]
        self.Z_syns = [np.array([1 + 4 * j + 2 * k, -1 + 4 * i + 2 * k])
                       for k in range(2) for j in range(2) for i in range(3)]
        # The part of the parity check circuit where noise can occur
        round_0 = {'RX': self.X_syns, 'RZ': self.Z_syns, 'IDLE_R': v in self.data}

        round_1 = {
            'CNOT': [(v, v + np.array([1, 1])) for v in self.X_syns if v[0] < 8] + [(v, v + np.array([1, 1])) for v in
                                                                                    self.Z_syns if v[1] < 8],
            'IDLE_CNOT': [np.array([0, 0]), np.array([8, 0]), np.array([0, 8])]}

        round_2 = {
            'CNOT': [(v, v + np.array([-1, 1])) for v in self.X_syns if v[0] > 0] + [(v, v + np.array([1, -1])) for v in
                                                                                     self.Z_syns if v[1] > 0],
            'IDLE_CNOT': [np.array([0, 0])]}

        round_3 = {
            'CNOT': [(v, v + np.array([1, -1])) for v in self.X_syns if v[0] < 8] + [(v, v + np.array([-1, 1])) for v in
                                                                                     self.Z_syns if v[1] < 8],
            'IDLE_CNOT': [np.array([8, 8])]}

        round_4 = {
            'CNOT': [(v, v + np.array([-1, -1])) for v in self.X_syns if v[0] > 0] + [(v, v + np.array([-1, -1])) for v in
                                                                                     self.Z_syns if v[1] > 0],
            'IDLE_CNOT': [np.array([8, 0]), np.array([0, 8]), np.array([8, 8])]}

        round_5 = {'MRX': self.X_syns, 'MRZ': self.Z_syns, 'IDLE_MR': self.data}

        # The part of the parity check circuit which is perfect
        round_6 = {
            'CNOT': [(v, v + np.array([1, 1])) for v in self.X_syns if v[0] < 8] + [(v, v + np.array([1, 1])) for v in
                                                                                    self.Z_syns if v[1] < 8]}

        round_7 = {
            'CNOT': [(v, v + np.array([-1, 1])) for v in self.X_syns if v[0] > 0] + [(v, v + np.array([1, -1])) for v in
                                                                                     self.Z_syns if v[1] > 0]}

        round_8 = {
            'CNOT': [(v, v + np.array([1, -1])) for v in self.X_syns if v[0] < 8] + [(v, v + np.array([-1, 1])) for v in
                                                                                     self.Z_syns if v[1] < 8]}

        round_9 = {
            'CNOT': [(v, v + np.array([-1, -1])) for v in self.X_syns if v[0] > 0] + [(v, v + np.array([-1, -1])) for v
                                                                                      in
                                                                                      self.Z_syns if v[1] > 0]}

        round_10 = {'RX': self.X_syns, 'RZ': self.Z_syns}

        # Store the circuit and fault locations
        self.circuit = [round_0, round_1, round_2, round_3, round_4, round_5, round_6, round_7, round_8, round_9, round_10]
        self.fault_locations = []
        for round_number in range(6):
            round = self.circuit[round_number]
            for fault_type in round.keys():
                target_qubits_list = round[fault_type]
                for qbs in target_qubits_list:
                    self.fault_locations.append((round_number, fault_type, qbs))

    def syn_diff(self, l, P):
        """Given an error (l, P), simulates the parity-check circuit followed by a perfect parity-check circuit,
        and returns the syndrome difference history as a (3, 4, 2) array plus a (4, 3, 2) array (for the
        X-syndromes and Z-syndromes, respectively)."""
        pass

    def error(self, l, P):
        """Given an error (l, P), returns the net effect on the data qubits as a shape (5, 5) array."""
        pass

    def edge_error_list(self, e1, e2):
        """Given vertices v1 = (x1, y1, t1) and v2 = (x2, y2, t2) with t1, t2 in {0, 1}, returns the set of
        faults (l, P) which trigger this edge."""
        pass
