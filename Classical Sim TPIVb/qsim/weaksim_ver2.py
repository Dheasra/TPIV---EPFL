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
        tmp = k.keys()
        # print('tmp =', tmp)
        for m in tmp:
            # print('m = ', m)
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
        if col[k] == 1 and (k != q): #Last condition: skip duplicates
            idx.append(k)
    outstr = ''
    for k in idx:
        outstr += str(k)
    return outstr

def entangled_qb(circ):
    #finds entangled qubits in the initial circuit
    out = []
    for k in range(len(circ[0])):
        if circ[0][k] == 'CX':
            out.append(circ[1][k])
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

def arr_is_zero(arr):
    for k in range(len(arr)):
        if arr[k] != 0:
            return False
    return True

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
            b[q] = (b[q]+1)%2 #TODO: This is the main issue right here
    
    #TODO: Implement these gates someday
    # elif gate == 'S':
    # elif gate == 'Z':
    # elif gate == 'Y':
    return R, T, c, b, d

def evol(circ, R, T, c, b, d):
    for g in range(len(circ[0])):
        R, T, c, b, d = upd_A(R, T, c, b, d, circ[0][g], circ[1][g])
    return R, T, c, b, d

def project(R,T,c,b,d, stb_st, qb):
    #Updates the subspace A according to the projectors in the decomposition of <T|
    #index of the ancilla qubit which will be projected
    qb = int(qb) 
    # print('qb' ,qb)
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
        #delete the columns and line in R where R[qb][l] = 1 
        idx = [l for l in range(len(R.mat[qb])) if R.mat[qb][l] == 1] #finding the index l of the 1 in R[qb]
        # print(idx)
        entqub = ''
        for k in idx:
            # print('T: ', T)
            #finding the entangled qubits
            tmp = entangled_idx(R.get_col(k), qb)
            entqub += tmp
        for l in range(len(entqub)):
            qrtl = int(entqub[l]) #control qubit of the CNOT linking the ancilla and the original qubit with T
            #updating A
            R.mat[qrtl] = add_vec(R.mat[qrtl], R.mat[qb], 2)
            # print('R: ', R.mat)
            T[qrtl] = (T[qrtl] + T[qb])%2
            # print('T: ', T)
            #updating q, l
            c.mat[qrtl] = add_vec(c.mat[qrtl], c.mat[qb], 2)
            b[qrtl] = (b[qrtl] + b[qb])%2
            d[qrtl] = (d[qrtl] + d[qb])%2

        #deleting the first column l where R[qb][l] = 1 
        if len(R.mat[0]) > 1 and len(idx) >0 and len(R.mat[0]) > (len(R.mat)-1):
            R.del_col(idx[0])
        del R.mat[qb]
        del T[qb]
        c.del_col(qb)
        del c.mat[qb]
        del b[qb]
        del d[qb]

    return R, T, c, b, d

def evol_proj(R, T, c, b, d, stb_states, qblist):
    if len(qblist) > 0:
        chi = len(stb_states) #approximate tensor rank
    else:
        chi = 1

    #creating the lists of R and T for each element of the stabilizer approximation of |A>
    Rlist = [deepcopy(R) for l in range(chi)]
    Tlist = [deepcopy(T) for l in range(chi)]
    clist = [deepcopy(c) for l in range(chi)]
    blist = [deepcopy(b) for l in range(chi)]
    dlist = [deepcopy(d) for l in range(chi)]

    #evolving the elements of R and T for each superposition 
    for l in range(len(qblist)):
        for k in range(len(stb_states)): 
            # print('k = ', k)
            # print('T[0] = ', Tlist[0])
            # print('T[1] = ', Tlist[1])
            Rlist[k], Tlist[k], clist[k], blist[k], dlist[k] = project(Rlist[k], Tlist[k], clist[k], blist[k], dlist[k], stb_states[k], qblist[l])
        for i in range(len(qblist)):
            qblist[i] = qblist[i]-1

    return Rlist, Tlist, clist, blist, dlist

def sample_shot(Rlist, Tlist, coefflist, ulist):
    N = len(Rlist[0].mat) #nbr of qubits in the final state
    chi = len(Rlist) #nbr of elements in the stabilizer decomposition of the magic state

    p = np.zeros(chi).tolist()
    for k in range(chi): #loop over stabilizer states in decomposition
        val = ''
        for l in range(N): #loop over each vector in the solution
            val += str(int(((mult_arr(Rlist[k].mat[l], ulist[k])) + Tlist[k][l])%2))
        p[k] = val
    return p

def sample_all(shots, Rlist, Tlist, coefflist):
    # chi = len(Rlist)
    nlist = [] #number of columns of each matrix in Rlist (it may not always be the same value for all matrices)
    for k in Rlist:
        nlist.append(len(k.mat[0]))
    
    # plist = np.zeros(shots).tolist() #list of bitstrings obtained by the circuit
    plist = [[],[]] #list of possible bitstrings obtained by the circuit 
    for k in range(shots):
        # generate uniformly random bitstrings u
        ulist = []
        for l in range(len(Rlist)):
            tmp = np.random.choice([0,1],nlist[l])
            ulist.append(tmp.tolist())
        #sample the probability distribution with the generated bitstrings
        temp = sample_shot(Rlist, Tlist, coefflist, ulist)
        for m in range(len(temp)):
            if not temp[m] in plist[m]:
                plist[m].append(temp[m])
    return plist
    # #reshape the bitstrings for each element in stabilizer decomp
    # nplist = [] #reshaped list of bitstrings
    # drac = [] #Count of the occurences of each bitstrings
    # for k in range(chi):
    #     tmp = plist.get_col(k)
    #     nplist.append(tmp) 
    #     #count the occurrences of each bitstring 
    #     addval = Counter(nplist[k])
    #     drac.append(addval)
    # #add dictionnaries together 
    # pfin = add_dict(drac, coefflist, True)
    # return pfin 

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

def compute_p(Rlist, Tlist, clist, blist, dlist, coefflist, shots): 
    n = len(Rlist[0].mat)
    # print('n = ', n)
    plist = [] #list of probabilities of outcome of each result for both Clifford circuits
    for k in range(n):
        plist.append([{'0': 1, '1': 1},{'0': 1, '1': 1}]) #The dictionaries correspond to each of the two circuits, each sublist correspond to one possible outcome
    kyes = ['0','1']
    # generate the list of all possible u bitstrings for every element of Rlist
    xlist = sample_all(shots, Rlist, Tlist, coefflist)
    # print('xlist = ', xlist)
    outcome = [{},{}]
    for k in range(len(Rlist)): #iterating on the stabilizer states
        xtmp = xlist[k]
        # print(xtmp)
        for p in kyes:
            for m in range(len(Rlist[k].mat)): #iterating on the qubits
                if not arr_is_zero(Rlist[k].mat[m]): #case qubit m is a superposition of |0> and |1>
                    #checking that the qubit is not entangled with a previous one
                    is_indep = True 
                    for r in range(m):
                        if mult_arr(Rlist[k].mat[m], Rlist[k].mat[r]) !=0:
                            is_indep = False
                    if is_indep:
                        plist[m][k][p] = plist[m][k][p]/np.sqrt(2)
                else: 
                    # print(xtmp[m])
                    if Tlist[k][m] == int(xtmp[m]): #useless for everything but clarity
                        plist[m][k][p] = plist[m][k][p]*1
                    else: 
                        plist[m][k][p] = 0
        #TODO: faire se multiplier les sous-listes de plist et s'additionner les dicitonnaires de plist
        #extracting the probability of each outcome from plist
        for a in range(pow(2,n)): #Iterating on all n bit outcomes
            binary = to_bin(a, n)
            globkey = ''
            prob = 1
            for m in range(n):
                cur_ky = str(binary[m])
                globkey += cur_ky 
                prob *= plist[m][k][cur_ky]
            if not globkey in xlist[k]: 
                prob = 0
            # print('globkey = ', globkey)
            phase = pow(-1, fct_q(clist[k], blist[k], binary))*pow(1j,fct_l(dlist[k], binary))
            # print('phase = ', phase)
            outcome[k].update({globkey: prob*phase})
            # print('outcome: ', k, outcome[k])
    # adding the probability amplitudes of every outcome together and taking the modulus
    real_out = add_dict(outcome, coefflist, True)
    # print('real_out =', real_out)
    return real_out

#Main code  
#===== Initial values
shots = 1000 #number of shots used for the determination of possible outcomes, should always be >2^Nqubits
#-----example_1
# Nqubits, circ = 1, [['H','T', 'X', 'H'], ['0', '0', '0', '0' ]]
#-----example_2 - simplest case
# Nqubits, circ = 1, [['H', 'T'], [ '0', '0']] 
# -----example_3 - 2 qubits full clifford 
# Nqubits, circ = 2, [['H', 'CX','X'], [ '0', '01','0']] 
#-----example_4 - reverse X-H from 1st example
# Nqubits, circ = 1, [['H','T', 'H', 'X'], ['0', '0', '0', '0' ]] 
#-----example_4, example with a controlled gate '01' means 0 is control qubit and 1 is target - Doesn't work
Nqubits, circ = 2, [['H', 'T', 'CX', 'H'], ['0', '0', '01', '0' ]]
#-----example_7, 2 independant qubits
# Nqubits, circ = 2, [['H','T', 'X', 'H', 'H','T', 'X', 'H'], ['0', '0', '0', '0', '1', '1', '1', '1']] 
#-----example_8, 4 independant qubits
# Nqubits, circ = 4, [['H','T', 'X', 'H', 'H','T', 'X', 'H', 'H','T', 'X', 'H', 'H','T', 'X', 'H'], ['0', '0', '0', '0', '1', '1', '1', '1','2', '2', '2', '2','3', '3', '3', '3' ]]

#===== Step 1: gadgetize
eq = entangled_qb(circ) #find which qubits are entangled together before gadgetization
circ, Nqubits, t, g_list = gagdetize(circ, Nqubits)

#===== step 2: Update A up to the projection onto the magic states 
R = Matrix(Nqubits, 1) #m=1 by default
T = np.zeros(Nqubits).tolist()
c = Matrix(Nqubits, Nqubits)
c.to_id()
b = np.ones(Nqubits).tolist()
d = np.zeros(Nqubits).tolist()
# print('R = ', R.mat)
R, T, c, b, d = evol(circ, R, T, c, b, d)

# print('R = ', R.mat)
# print('T = ', T)
# print('c = ', c.mat)
# print('b = ', b)
# print('d = ', d)

#===== step 3: Approximate magic states as linear combination of stabilizer states
#states
if len(g_list)>0:
    stab= [[''],['X']]
    coeff = [1/np.sqrt(2), cmath.exp(math.pi*1j/4)/np.sqrt(2)]
else: 
    stab = [['']]
    coeff = [1]
# stab= [['H'],['X']] #Each "row" (i.e. sub-list) represents a stabilizer state present in the sum
                    #Each stabilizer is represented by the clifford circuit generating it from |0>, with the gate at index 0 being the first to be applied, 1 the second, etc


#===== step 4: Update R, T with the projectors onto each element of the linear combination
Rlist, Tlist, clist, blist, dlist = evol_proj(R, T, c, b, d, stab, g_list)

# print('Rlist[0]: ', Rlist[0].mat) #note to self: in example_1, notice that this is a 1x2 matrix, but it doesn't really matter for the first element is 0 and thus this matrix represents the state |0>+|1> as it should
# print('Tlist[0]: ', Tlist[0])
# print('clist[0]: ', clist[0].mat) 
# print('blist[0]: ', blist[0])
# print('dlist[0]: ', dlist[0])
# print('Rlist[1]: ', Rlist[1].mat) 
# print('Tlist[1]: ', Tlist[1])
# print('clist[1]: ', clist[1].mat) 
# print('blist[1]: ', blist[1])
# print('clist[1]: ', dlist[1])

#===== step 5: sample p for each elements of the linear combination
psample = sample_all(shots, Rlist, Tlist, coeff) 
# print(psample)
out = compute_p(Rlist, Tlist, clist, blist, dlist, coeff, shots)
print(out)

#Visualization of the result
plt.bar(range(len(out)), list(out.values()), align='center')
plt.xticks(range(len(out)), list(out.keys()))
# plt.hist(p)
plt.show()

# print('circuit: ', circ)
# print('Nq: ', Nqubits)
# print('T count: ', t)
# print('g_list: ', g_list)


# tmp1  = [1, 1, 0, 0]
# str1 = entangled_idx(tmp1)
# print(str1)