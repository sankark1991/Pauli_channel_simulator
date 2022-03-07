

class BaseParityCheck:
    """Rotated surface code with distance 5."""
    def __init__(self):
        """Gets and stores the full list of error locations/types (l, P).
        Pre-computes and stores syn_diff(l, P) and error(l, P) for each one.
        Then for each edge (v1, v2), stores the list of faults (l, P) which trigger that edge."""
        pass

    def circuit(self):
        """Encodes the parity check circuit, and also stores the list of error locations/types (l, P)."""
        pass

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
