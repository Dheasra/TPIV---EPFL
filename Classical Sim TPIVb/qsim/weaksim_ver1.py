import math
import cmath
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib
import random
from copy import deepcopy
from collections import Counter

#Classes
class Matrix:
    mat : list
    def __init__(self, N = 1,m = 1):
        self.mat = np.zeros((N,m))
        self.mat = self.mat.tolist()

    def append_col(self, addcol):
        for k in range(len(self.mat)):
            self.mat[k].append(addcol[k])
    
    def is_zero(self, idx = -1):
        if idx == -1: #check all matrix
            for k in range(len(self.mat)):
                for l in range(len(self.mat[k])):
                    if self.mat[k][l] != 0:
                        return False
        elif idx > -1: #check only column with index = idx
            for k in range(len(self.mat)):
                if self.mat[k][idx] != 0:
                    return False
        return True

    def del_col(self, idx):
        for k in range(len(self.mat)):
            del self.mat[k][idx]
    
    def get_col(self, idx):
        col = []
        for k in range(len(self.mat)):
            col.append(self.mat[k][idx])
        return col

class bin_arr:
    #Binary array
    arr : list
    def __init__(self, new_arr):
        # self.arr = np.zeros(len(new_arr)).tolist()
        # for k in range(len(new_arr)):
        #     self.arr[k] = new_arr[k]
        self.arr = []
        for k in new_arr:
            self.arr.append(k)
    
    @classmethod
    def init_zero_arr(self, N):
        return bin_arr(np.zeros(N).tolist())

    def bin_add(self, add_arr):
        for k in range(len(self.arr)):
            self.arr[k] = (self.arr[k] + add_arr[k])%2

#Utilities
def mult_arr(arr1, arr2):
    #Scalar product of two arrays (lists)
    rslt = 0
    for k in range(len(arr1)):
        rslt += arr1[k]*arr2[k]
    return rslt

def normalize(arr):
    #Normalize an array
    norm = 0
    for k in arr:
        norm += pow(abs(k), 2)
    norm = pow(norm, 1/2)
    if norm == 0:
        return arr
    else: 
        return [k/norm for k in arr]

def add_vec(arr1, arr2, mod = 0):
    #simple addition of 2 arrays (lists)
    arr_f = []
    for k in range(len(arr1)):
        arr_f.append(arr1[k] + arr2[k])
        if mod != 0:
            arr_f[k] = arr_f[k] % mod

    return arr_f

def add_dict(dlist, coeffs = 0, modulus = False):
    #Adds the values of 2 dictionaries according to their keys
    n = len(dlist)
    if coeffs == 0:
        coeffs = np.ones(n).tolist()
    key = []
    for k in dlist: #get the keys of all the dictionaries
        key.append(k.keys())
    key = flatten(key) 
    sumd = {} #sum of the dictionnaries
    for k in key:
        sumd.update({k : 0})
        for l in range(len(coeffs)):
            sumd[k] += coeffs[l]*dlist[l][k]
    if modulus:
        for k in key:
            sumd[k] = pow(abs(sumd[k]),2)
    return sumd

def entangled_idx(col, q):
    #Finds the indices of qubits entangled with qubit q and return them. If there are none returns empty string
    idx = []
    for k in range(len(col)):
        if col[k] == 1 and (k != q): #Last condition: skip duplicates
            idx.append(k)
    outstr = ''
    for k in range(len(idx)):
        outstr += str(k)
    return outstr

def flatten(mat):
    #Flattens the matrix row wise, so that the end matrix has only 1 row
    out = []
    for k in range(len(mat)):
        for l in mat[k]:
            out.append(l)
    return out

#Core computations
def gagdetize(circuit, Nq):
    #gadgetizes the circuit 
    t = 0 #T count
    g_list = []#list of qubits added by gadgetization
    #Changing the T gates to the T Gadget
    for k in range(len(circuit[0])):
        if circuit[0][k] == 'T':
            t += 1 
            Nq += 1 #adding a qubit for every T gate
            g_list.append(Nq-1) #adding the newly added qubit to the list

            circuit[0][k] = 'CX' #Replacing the gate 
            control = int(circuit[1][k]) #the control qubit is where the T gate was
            circuit[1][k] = str(control) + str(Nq-1) #changing the qubits where the gate is applied
    return circuit, Nq, t, g_list

def upd_A(R, T, gate, qb):
    #Updates the subspace A according to the gate "gate" applied on qubit qb
    if gate =='H':
        new_col = True 
        q = int(qb)
        #Verifying that R must have a supplementary column  and verifying if the qubit is entangled with another
        addcol = bin_arr.init_zero_arr(len(R.mat)) #column to be added to R
        addcol.arr[q] = 1
        for l in range(len(R.mat[q])):
            if R.is_zero(l): #Checking that column l contains only zeros (if it does => no supplementary column needed)
                new_col = False
                no_col = l
            if R.mat[q][l] == 1: #Checking if qubit q is entangled with any other
                addcol.bin_add(deepcopy(R.get_col(l))) #adding the columns corresponding to the entanglement
        if new_col:
            R.append_col(deepcopy(addcol.arr))
        else:
            R.mat[q][no_col] = 1
        
    elif gate == 'CX':
        c = int(qb[0])
        t = int(qb[1])
        for k in range(len(R.mat[c])):
            R.mat[t][k] = (R.mat[t][k] +  R.mat[c][k])%2
        T[t] = (T[t] + T[c])%2

    elif gate == 'X':
        q = int(qb)
        T[q] = (T[q]+1)%2
                    
    # elif gate == 'S':
    # elif gate == 'Z':
    # elif gate == 'Y':
    return R, T

def evol(circ, R, T):
    for g in range(len(circ[0])):
        R, T = upd_A(R, T, circ[0][g], circ[1][g])
    return R, T

def project(R,T, stb_st, qb):
    #Updates the subspace A according to the projectors in the decomposition of <T|
    qb = int(qb) 
    #check if the projection gives zero probability
    p_zero = False
    if R.is_zero() and not ('H' in stb_st): #both R and R_stab are zero
        for k in stb_st:
            if k == 'X' and T[qb] == 0: #if k=X then T_stab=1 , and if T_stab != T[qb] then the projection gives amplitude 0
                p_zero = True
            elif not (k=='X') and T[qb] == 1: #this is the case where T_stab = 0 and and T[qb] = 1
                p_zero = True
    if p_zero: 
        R = Matrix()
        R.mat[0][0] = math.nan 
        T = np.zeros(1).tolist()
        T[0] = math.nan
    else: 
        for k in stb_st: #updating the state with the Clifford circuit of the stabilizer state before projection
            R, T = upd_A(R, T, k, qb)
        #delete the columns and line in R where R[qb][l] = 1 
        idx = [l for l in range(len(R.mat[qb])) if R.mat[qb][l] == 1] #finding the index l of the 1 in R[qb]
        entqub = ''
        for k in idx:
            #finding the entangled qubits
            tmp = entangled_idx(R.get_col(k), qb)
            entqub += tmp
        for l in range(len(entqub)):
            qrtl = int(entqub[l]) #control qubit
            R.mat[qrtl] = add_vec(R.mat[qrtl], R.mat[qb], 2)
            T[qrtl] = (T[qrtl] + T[qb])%2
        #deleting the first column l where R[qb][l] = 1 
        if len(R.mat[0]) > 1 and len(idx) >0:
            R.del_col(idx[0])
        del R.mat[qb]
        del T[qb]

    return R, T 

def evol_proj(R, T, stb_states, qblist):
    chi = len(stb_states) #approximate tensor rank

    #creating the lists of R and T for each element of the stabilizer approximation of |A>
    Rlist = [deepcopy(R) for l in range(chi)]
    Tlist = [deepcopy(T) for l in range(chi)]

    #evolving the elements of R and T for each superposition 
    for l in range(len(qblist)):
        for k in range(len(stb_states)): 
            # print('k = ', k)
            # print('T[0] = ', Tlist[0])
            # print('T[1] = ', Tlist[1])
            Rlist[k], Tlist[k] = project(Rlist[k], Tlist[k], stb_states[k], qblist[l])

    return Rlist, Tlist 

def sample_shot(Rlist, Tlist, clist, ulist): #TODO: corriger, le problème vient du val += ... probablement -> normalement c'est corrigé
    N = len(Rlist[0].mat) #nbr of qubits in the final state
    chi = len(Rlist) #nbr of elements in the stabilizer decomposition of the magic state

    p = np.zeros(chi).tolist()
    for k in range(chi): #loop over stabilizer states in decomposition
        val = ''
        for l in range(N): #loop over each vector in the solution
            val += str(int(((mult_arr(Rlist[k].mat[l], ulist[k])) + Tlist[k][l])%2))
        p[k] = val
    return p

def sample_all(shots, Rlist, Tlist, clist):
    chi = len(Rlist)

    nlist = [] #number of columns of each matrix in Rlist (it may not always be the same value for all matrices)
    for k in Rlist:
        nlist.append(len(k.mat[0]))
    
    # plist = np.zeros(shots).tolist() #list of bitstrings obtained by the circuit
    plist = Matrix(shots) #list of bitstrings obtained by the circuit
    for k in range(shots):
        # generate uniformly random bitstrings u
        ulist = []
        for l in range(len(Rlist)):
            tmp = np.random.choice([0,1],nlist[l])
            ulist.append(tmp.tolist())
        #sample the probability distribution with the generated bitstrings
        plist.mat[k] = sample_shot(Rlist, Tlist, clist, ulist)
    #reshape the bitstrings for each element in stabilizer decomp
    nplist = [] #reshaped list of bitstrings
    drac = [] #Count of the occurences of each bitstrings
    for k in range(chi):
        tmp = plist.get_col(k)
        nplist.append(tmp) 
        #count the occurrences of each bitstring 
        addval = Counter(nplist[k])
        drac.append(addval)
    #add dictionnaries together 
    pfin = add_dict(drac, clist,False)
    return pfin 

#Main code  
#===== Initial values
shots = 1024 #number of shots for the sampling
Nqubits, circ = 1, [['H','T', 'X', 'H'], ['0', '0', '0', '0' ]] #example_1
# Nqubits, circ = 1, [['H', 'T'], [ '0', '0']] #example_2 - simplest case
#Nqubits, circ = 2, [['H', 'CX','X'], [ '0', '01','0']] #example_3 - full clifford
#Nqubits, circ = 6, [['H', 'CX','CX', 'CX','H', 'CX'], [ '0', '01','02', '23','4', '45']] #example_3bis - full clifford
#Nqubits, circ = 2, [['H', 'T', 'CX', 'H'], ['0', '0', '01', '1' ]] #example_4, example of a controlled gate '01' means 0 is control qubit and 1 is target
# Nqubits, circ = 2, [['H', 'CX','T'], [ '0', '01','0']] #example_5 - NF style circuit

#===== Step 1: gadgetize
circ, Nqubits, t, g_list = gagdetize(circ, Nqubits)

#===== step 2: Update A up to the projection onto the magic states 
R = Matrix(Nqubits, 1) #m=1 by default
T = np.zeros(Nqubits).tolist()
# print('R = ', R.mat)
R, T = evol(circ, R, T)

print('R = ', R.mat)
print('T = ', T)

#===== step 3: Approximate magic states as linear combination of stabilizer states
#states
# stab= [[''],['X']]
# coeff = [1, cmath.exp(math.pi*1j/4)]
stab= [[''],['X']]
coeff = [1/2, 1/2]
# stab = [['']]
# coeff = [1]
# stab= [['H'],['X']] #Each "row" (i.e. sub-list) represents a stabilizer state present in the sum
                    #Each stabilizer is represented by the clifford circuit generating it from |0>, with the gate at index 0 being the first to be applied, 1 the second, etc


#===== step 4: Update R, T with the projectors onto each element of the linear combination
Rlist, Tlist = evol_proj(R, T, stab, g_list)

# print('Rlist[0]: ', Rlist[0].mat) #note to self: in example_1, notice that this is a 1x2 matrix, but it doesn't really matter for the first element is 0 and thus this matrix represents the state |0>+|1> as it should
# print('Tlist[0]: ', Tlist[0])
# print('Rlist[1]: ', Rlist[1].mat) 
# print('Tlist[1]: ', Tlist[1])

#===== step 5: sample p for each elements of the linear combination
p = sample_all(shots, Rlist, Tlist, coeff) #TODO: corriger - cause probable: on néglige totalement les phases des états individuels qui peuvent s'annuler et s'additionner -> voir le dernier papier de bravyi 
print(p)

#Visualization of the result
plt.bar(range(len(p)), list(p.values()), align='center')
plt.xticks(range(len(p)), list(p.keys()))
# plt.hist(p)
plt.show()

# print('circuit: ', circ)
# print('Nq: ', Nqubits)
# print('T count: ', t)
# print('g_list: ', g_list)


# tmp1  = [1, 1, 0, 0]
# str1 = entangled_idx(tmp1)
# print(str1)