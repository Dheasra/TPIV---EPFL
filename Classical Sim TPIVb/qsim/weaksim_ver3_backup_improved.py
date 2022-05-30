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

    def set_col(self, col, idx):
        for k in range(len(self.mat)):
            self.mat[k][idx] = col[k]
    
    def to_id(self):
        if len(self.mat) == len(self.mat[0]):
            for k in range(len(self.mat)): #Purging the matrix
                for l in range(len(self.mat[k])):
                    self.mat[k][l]  = 0
                self.mat[k][k] = 1 #Replacing the main diagonal with 1
        else: 
            print('Warning: Cannot transform non square matrix to identity')

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

def normalize(dico):
    #Normalize a dictionary
    total = 0
    #Find the total of the elements in the dictionary
    for k in dico.keys():
        total += dico[k]
    #Now divide each element in the dictionary by the total
    for k in dico.keys():
        dico[k] /= total
    return dico
    
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
        tmp = k.keys()
        for m in tmp:
            if not m in key:
                key.append(m)
    sumd = {} #sum of the dictionnaries
    for k in key:
        sumd.update({k : 0})
        for l in range(len(coeffs)):
            if k in dlist[l]:
                sumd[k] += coeffs[l]*dlist[l][k]
    if modulus:
        for k in key:
            sumd[k] = pow(abs(sumd[k]),2)
    return sumd

def entangled_idx(col, q):
    #Finds the indices of qubits entangled with qubit q and return them. If there are none returns empty string
    idx = []
    for k in range(len(col)):
        if col[k] == 1 and (k != q) and not (k in idx): #Last condition: skip duplicates
            idx.append(k)
    outstr = ''
    for k in idx:
        outstr += str(k)
    return outstr

def entangled_qb(circ, Nq):
    #finds entangled qubits in the initial circuit
    entqub = []
    for k in range(len(circ[0])):
        if circ[0][k] == 'CX':
            entqub.append(circ[1][k])
    out = entqub
    # out = []
    # #group all entangled qubits together
    # for i in range(Nq):
    #     for k in entqub:
    #         if str(i) in k:

    return out

def to_bin(n, bits = 0):
    b = []
    while n>0:
        d = n%2
        b.append(d)
        n=n//2
    if len(b) < bits:
        for a in range(bits-len(b)):
            b.append(0)
    b.reverse()
    return b

def str2list(in_str):
    out = []
    for k in in_str:
        out.append(int(k))
    return out

def arr_is_zero(arr):
    for k in range(len(arr)):
        if arr[k] != 0:
            return False
    return True

def count_occurrences(arr, cmat, bvec, dvec):
    #count the occurences of elements in arr with phase according to c,b,d (functions q and l)
    dico = {}
    for k in arr: #TODO: finir
        if k in dico.keys():
            dico[k] += pow(-1,fct_q(cmat, bvec, str2list(k)))*pow(1j, fct_l(dvec, str2list(k)))
        else: 
            dico.update({k: pow(-1,fct_q(cmat, bvec, str2list(k)))*pow(1j, fct_l(dvec, str2list(k)))})
    return dico

def clear_duplicates(arr):
    #removes all duplicate elements in string arr
    clr = []
    for k in arr:
        if not k in clr:
            clr.append(k)
    #transforms clr into a string
    out = ''
    for k in clr:
        out += k
    return out

def to_str(lst):
    #converts a list to a string without blanks in between elements
    out = ''
    for k in lst:
        out += str(k)
    return out

def estim_error(dicolist):
    #Averages the values in dicolist and gives the standard error
    key = []
    for k in dicolist: #get the keys of all the dictionaries
        tmp = k.keys()
        for m in tmp:
            if not m in key:
                key.append(m)
    out = {} #averaged results
    err = {} #standard error
    for k in key:
        val = []
        for l in dicolist:
            val.append(l[k])
        #averaging: 
        out.update({k: np.mean(val)})
        err.update({k: np.std(val)})
    return out, err


#Core computations
def gagdetize(circuit, Nq):
    #gadgetizes the circuit 
    t = 0 #T count
    g_list = [] #list of qubits added by gadgetization
    ctrl_list = [] #list of entangled qubits with the new ancilla qubits
    #Changing the T gates to the T Gadget
    for k in range(len(circuit[0])):
        if circuit[0][k] == 'T':
            ctrl_list.append(int(circuit[1][k]))
            t += 1 
            Nq += 1 #adding a qubit for every T gate
            g_list.append(Nq-1) #adding the newly added qubit to the list

            circuit[0][k] = 'CX' #Replacing the gate 
            control = int(circuit[1][k]) #the control qubit is where the T gate was
            circuit[1][k] = str(control) + str(Nq-1) #changing the qubits where the gate is applied
    return circuit, Nq, t, g_list, ctrl_list

def upd_A(R, T, c, b, d, gate, qb, xstab = False):
    #Updates the triplet (A,q,l) according to the gate "gate" applied on qubit qb
    #Recall that A is determined by R and T, q by c and b and l by d
    if gate =='H':
        #evolving A
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
        #evolving c
        #Constructing the vector v to add to row q of c
        v = np.zeros(len(R.mat)).tolist() #v is a vector containing 1 at every index corresponding to qubits entangled to q
        for k in range(len(R.mat[q])):
            if R.mat[q][k] == 1:
                tmp_col = R.get_col(k)
                for l in range(len(tmp_col)):
                    if tmp_col[l] == 1:
                        v[l] = 1
        #adding it to row q of c
        c.set_col(add_vec(deepcopy(c.mat[q]), v, 2), q)
        #evolving b
        b[q] = (b[q] + T[q])%2
        
    elif gate == 'CX':
        ctrl = int(qb[0])
        trgt = int(qb[1])
        for k in range(len(R.mat[ctrl])):
            R.mat[trgt][k] = (R.mat[trgt][k] +  R.mat[ctrl][k])%2
        T[trgt] = (T[trgt] + T[ctrl])%2

    elif gate == 'X':
        q = int(qb)
        T[q] = (T[q]+1)%2
        # if not arr_is_zero(R.mat[q]): #TODO:doesn't work
        #     b[q] = (b[q]+1)%2 
        if xstab: #TODO: replace this with a system that works in all cases were X follows a H, not only for the stabilzer
            b[q] = (b[q]+1)%2 #TODO: This might be an issue
    
    #TODO: Implement these gates someday
    # elif gate == 'S':
    # elif gate == 'Z':
    # elif gate == 'Y':
    return R, T, c, b, d

def evol(circ, R, T, c, b, d):
    for g in range(len(circ[0])):
        R, T, c, b, d = upd_A(R, T, c, b, d, circ[0][g], circ[1][g])
    return R, T, c, b, d

def project(R,T,c,b,d, stb_st, qb, ctrl):
    #Updates the subspace A according to the projectors in the decomposition of <T|
    #index of the ancilla qubit which will be projected
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
            R, T, c, b, d = upd_A(R, T, c, b, d, k, qb, True)
        idx = [l for l in range(len(R.mat[qb])) if R.mat[qb][l] == 1] #finding the index l of the 1 in R[qb]
        # print('idx =',idx, '|| qb =', qb)
        entqub = ''
        for k in idx:
            # print('T: ', T)
            #finding the entangled qubits
            tmp = entangled_idx(R.get_col(k), qb)
            entqub += tmp
        # print('entqub before ', entqub)
        entqub = clear_duplicates(entqub) #TODO: corriger 
        # print('entqub after ', entqub)
        for l in range(len(entqub)):
            qrtl = int(entqub[l]) #control qubit of the CNOT linking the ancilla and the original qubit with T
            #updating A
            R.mat[qrtl] = add_vec(R.mat[qrtl], R.mat[qb], 2)
            T[qrtl] = (T[qrtl] + T[qb])%2
            # print( 'l: ',l,'|| qrtl = ', qrtl, '|| ctrl = ', ctrl)
            #updating q, l
            c.mat[qrtl] = add_vec(c.mat[qrtl], c.mat[qb], 2)
            b[qrtl] = (b[qrtl] + b[qb])%2
            d[qrtl] = (d[qrtl] + d[qb])%2

        #deleting the first column l where R[qb][l] = 1 
        if len(R.mat[0]) > 1 and len(idx) >0:
            R.del_col(idx[0])
        del R.mat[qb]
        del T[qb]
        c.del_col(qb)
        del c.mat[qb]
        del b[qb]
        del d[qb]

    return R, T, c, b, d

def evol_proj(R, T, c, b, d, stb_states, qblist, ctrllist, t, coefflist):
    chi = len(stb_states) #approximate tensor rank
    nbr_of_list = pow(chi,t) #number of lists

    #creating the lists of R and T for each element of the stabilizer approximation of |A>
    Rlist = [deepcopy(R) for l in range(nbr_of_list)]
    Tlist = [deepcopy(T) for l in range(nbr_of_list)]
    clist = [deepcopy(c) for l in range(nbr_of_list)]
    blist = [deepcopy(b) for l in range(nbr_of_list)]
    dlist = [deepcopy(d) for l in range(nbr_of_list)]

    #creating the list of which stabilizer projectors are applied to which elements of the lists above #TODO: generalize it for any tensor rank, not just 2
    choice_list = [to_str(to_bin(l, t)) for l in range(nbr_of_list)]
    #computing the coefficient associated with the 
    newcoefflist = []
    for k in choice_list:
        newcoeff = 1
        for l in k:
            newcoeff *= coefflist[int(l)]
        newcoefflist.append(newcoeff)

    #evolving the elements of R and T for each superposition 
    for l in range(len(qblist)):
        # print('l =', l)
        for m in range(nbr_of_list): 
            k = int(choice_list[m][l])
            # print('m = ', m,'k = ', k)
            Rlist[m], Tlist[m], clist[m], blist[m], dlist[m] = project(Rlist[m], Tlist[m], clist[m], blist[m], dlist[m], stb_states[k], qblist[l], ctrllist[l])
        for i in range(len(qblist)):
            qblist[i] = qblist[i]-1
    return Rlist, Tlist, clist, blist, dlist, newcoefflist

def sample_shot(Rlist, Tlist, ulist): 
    N = len(Rlist[0].mat) #nbr of qubits in the final state
    chi = len(Rlist) #nbr of elements in the stabilizer decomposition of the magic state

    p = np.zeros(chi).tolist()
    for k in range(chi): #loop over stabilizer states in decomposition
        val = ''
        for l in range(N): #loop over each vector in the solution
            val += str(int(((mult_arr(Rlist[k].mat[l], ulist[k])) + Tlist[k][l])%2))
        p[k] = val
    return p

def sample_all(shots, Rlist, Tlist, clist, blist, dlist, coefflist):
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
        plist.mat[k] = sample_shot(Rlist, Tlist, ulist)
    #reshape the bitstrings for each element in stabilizer decomp
    nplist = [] #reshaped list of bitstrings
    drac = [] #Count of the occurences of each bitstrings
    for k in range(chi):
        tmp = plist.get_col(k)
        nplist.append(tmp) 
        #count the occurrences of each bitstring 
        addval = count_occurrences(nplist[k], clist[k], blist[k], dlist[k])
        drac.append(addval)
    #add dictionnaries together 
    pfin = add_dict(drac, coefflist,True)
    return pfin

def fct_q(c_mat, b_vec, bits):
    #quadratic function q used to determine the phase
    #multiplication x^t*C*x 
    tmp = []
    for i in range(len(c_mat.mat)): 
        tmp.append(mult_arr(c_mat.mat[i], bits))
    xcx = (mult_arr(tmp, bits))%2
    #multiplication b^t*x
    bx = (mult_arr(b_vec, bits))%2
    return (xcx+bx)%2

def fct_l(d_vec, bits):
    #linear function q used to determine the individual phases
    return (mult_arr(d_vec, bits))%2

#Main code  
#===== Initial values - Note: Circuits do not support S,Z and Y gates (yet)
shots = 1024 #number of shots used for the determination of possible outcomes, should always be >2^Nqubits
repet = 16 #number of repetitions of the sampling method to estimate the error on the results
#-----example_1
# Nqubits, circ = 1, [['H','T', 'X', 'H'], ['0', '0', '0', '0' ]]
#-----example_2 - simplest case
# Nqubits, circ = 1, [['H', 'T'], [ '0', '0']] 
# -----example_3 - 2 qubits full clifford 
# Nqubits, circ = 2, [['H', 'CX','X'], [ '0', '01','0']] 
# -----example_4, example of a controlled gate '01' means 0 is control qubit and 1 is target
# Nqubits, circ = 2, [['H', 'T', 'CX', 'H'], ['0', '0', '01', '0' ]] 
#-----example_5 - reverse X-H from 1st example
# Nqubits, circ = 1, [['H','T', 'H', 'X'], ['0', '0', '0', '0' ]] 
#-----example_7, 2 independant qubits
# Nqubits, circ = 2, [['H','T', 'X', 'H', 'H','T', 'X', 'H'], ['0', '0', '0', '0', '1', '1', '1', '1']] 
#-----example_8, 4 independant qubits
Nqubits, circ = 4, [['H','T', 'X', 'H', 'H','T', 'X', 'H', 'H','T', 'X', 'H', 'H','T', 'X', 'H'], ['0', '0', '0', '0', '1', '1', '1', '1','2', '2', '2', '2','3', '3', '3', '3' ]]
#-----example_8, 2 entangled qubits with multiple T gates
# Nqubits, circ = 2, [['H','T','H', 'CX','T', 'X'], ['0', '0', '0','01', '1', '1']] #TODO: faire une gestion des NaN dans le sampling

#===== Step 1: gadgetize
eq = entangled_qb(circ, Nqubits)
circ, Nqubits, t, g_list, ctrl_list = gagdetize(circ, Nqubits)

#===== step 2: Update A up to the projection onto the magic states 
R = Matrix(Nqubits, 1) #m=1 by default
T = np.zeros(Nqubits).tolist()
c = Matrix(Nqubits, Nqubits)
c.to_id()
b = np.ones(Nqubits).tolist()
d = np.zeros(Nqubits).tolist()
R, T, c, b, d = evol(circ, R, T, c, b, d)

# print('R = ', R.mat)
# print('T = ', T)
# print('c = ', c.mat)
# print('b = ', b)
# print('d = ', d)

#===== step 3: Approximate magic states as linear combination of stabilizer states
#states
stab= [[''],['X']]
# coeff = [1, cmath.exp(math.pi*1j/4)] #unnormalized coefficients (debugging)
coeff = [1/np.sqrt(2), cmath.exp(math.pi*1j/4)/np.sqrt(2)] #normalized coefficients

#===== step 4: Update R, T with the projectors onto each element of the linear combination
Rlist, Tlist, clist, blist, dlist, newcoeff = evol_proj(R, T, c, b, d, stab, g_list, ctrl_list, t, coeff)
# print('newcoeff', newcoeff)
# for k in range(len(Rlist)):
#     print('Rlist[',k,']: ',Rlist[k].mat)
#     print('Tlist[',k,']: ',Tlist[k])

#===== step 5: sample p for each elements of the linear combination
tmp_out= [{} for k in range(repet)]
for k in range(repet): #redo the sampling multiple times to have an estimate of the error
    tmp_out[k] = sample_all(shots, Rlist, Tlist, clist, blist, dlist, newcoeff) 
    tmp_out[k] = normalize(tmp_out[k])
out, err = estim_error(tmp_out)
print(out)
print(err)

#Visualization of the result
plt.bar(range(len(out)), list(out.values()),yerr= list(err.values()),align='center')
plt.xticks(range(len(out)), list(out.keys()))
# plt.hist(p)
plt.show()
