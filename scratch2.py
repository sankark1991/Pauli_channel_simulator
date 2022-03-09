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
from surface_code.BaseParityCheck import BaseParityCheck
from surface_code.ParityCheck import ParityCheck

load_data = False
if load_data:
    f = open("base_pc_data.pckl", 'rb')
    small_code = pickle.load(f)
    f.close()
else:
    small_code = BaseParityCheck()
    f = open("base_pc_data.pckl", 'wb')
    pickle.dump(small_code, f)
    f.close()


### Function for taking an RQP for the big code, and going down to an RQP for the small code
code_dist = 31
N_rounds = 5
PC = ParityCheck(code_dist, N_rounds, small_code)
error, sdiff = PC.sample(p=0.1)
print(error)
print(sdiff)
raise InterruptedError
