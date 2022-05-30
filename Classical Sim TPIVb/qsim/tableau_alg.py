import math
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib

#TODO:many things

def init_state(Nq, t = 0): #If t is given and bigger than 0: forward T gagdets are used, otherwise: reverse T gagdets are used  %TODO: implémenter le forward
    ste = np.zeros((2*Nqubits, 2*Nqubits + 1))
    for k in range(2*Nqubits):
        ste[k][k] = 1
    return ste

def upd_state(state, gate,qb):#TODO: finir
    Nq = len(state)

    if gate == 'H':
        a = int(qb)
        for k in range(Nq):
            r = state[k][2*Nq]
            x = state[k][a-1]
            z = state[k][a-1+Nq/2]

            r = (r + x*z)%2 
            tmp = x #x_k,a
            x = z 
            z = tmp
        return state
    elif gate == 'S':
        a = int(qb)
        for k in range(Nq):
            r = state[k][2*Nq]
            x = state[k][a-1]
            z = state[k][a-1+Nq/2]

            r = (r + x*z)%2 
            z = (x+z)%2
        return state
    elif gate == 'CX':
        a = int(qb[0]) #control
        b = int(qb[1]) #target
        for k in range(Nq):
            r = state[k][2*Nq]
            xa = state[k][a-1]
            za = state[k][a-1+Nq/2]
            xb = state[k][b-1]
            zb = state[k][b-1+Nq/2]

            r = (r + xa*zb*(xb+za+1))%2 
            xb = (xb+xa)%2
            za = (za+zb)%2
    elif gate == 'X':

    elif gate == 'Z':

    elif gate == 'Y':


    

#=== Initial values
Nqubits = 2 
circ = [['H', 'T', 'CX', 'H'], ['0', '0', '01', '1' ]] #circuit, example of a controlled gate '01' means 0 is control qubit and 1 is target

#Step 1: gadgetize
circ, Nqubits, t = gagdetize(circ, Nqubits)
print('circuit: ', circ)
print('Nq: ', Nqubits)
print('T count: ', t)

#step 2: Update A   
R = Matrix(Nqubits, 1) #m=1 by default
T = np.array(Nqubits)