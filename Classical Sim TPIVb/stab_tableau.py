import math
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib

#nbr of qubits
N = 10

class circuit:

    def __init__(self, N):
        #init tableau to initial state |0^N>
        self.tableau = np.zeros((2*N, 2*N+1))
        for j in range(N):
            self.tableau[j][j] = 1

    def CNOT(self, ctrl, trgt):
        for k in range(length(self.tableau)/2-1):
            




