def project(R,T, qb, Rproj, Tproj):
    #check if the projection gives zero probability
    p_zero = True
    if R.is_zero() and Rproj.is_zero():
        for k in Tproj:
            for l in T:
                if k==l:
                    p_zero = False 
    if p_zero:
        R = Matrix()
        R[0][0] = math.nan 
        T = np.zeros(1).tolist()
        T[0] = math.nan
    else: 
        qb = int(qb)
        del R.mat[qb] 
        del T[qb]
        # if len(R.mat[0]) > 1: 
        # else:
#============================================================

#first element 
Rp1 = Matrix() 
Tp1 = np.zeros(Nqubits).tolist()
Rp1.mat[0][0] = 1

#second element
Rp2 = Matrix() 
Tp2 = np.zeros(Nqubits).tolist()
Tp2[0] = 1

#============================================================

def is_zero_list(lst):
    for k in lst:
        if k != 0:
            return False
    return True

#============================================================

#backup 5/5/21 before I update the behaviour of H on R,T 
def upd_A(R, T, gate, qb):#TODO: checker que mon idée marche sur des exemples, puis l'implémenter si oui /maybe add list_qb is the list of qubits on which a hadamard has been applied
    if gate =='H':
        new_col = False #Checks if the qubit is entangled with another, or if it already had an hadamard applied
        q = int(qb)
        for l in range(len(R.mat)):
            for k in R.mat[l]:
                if k == 1:
                    new_col = True
        if new_col:
            addcol = np.zeros(len(R.mat))
            addcol[int(q)] = 1
            R.append_col(addcol)
        else:
            R.mat[q][0] = 1

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

#============================================================

    def gen_rdm_bitstr(n):
    bstr = []
    for i in range(n):
        bstr.append(random.randint(0,1))
    return bstr

#============================================================

def sample_all(shots, Rlist, Tlist, clist):
    N = len(Rlist[0].mat) #nbr of qubits

    print('Rlist[0]: ', Rlist[0].mat) #note to self: in example_1, notice that this is a 1x2 matrix, but it doesn't really matter for the first element is 0 and thus this matrix represents the state |0>+|1> as it should
    print('Tlist[0]: ', Tlist[0])
    print('Rlist[1]: ', Rlist[1].mat) 
    print('Tlist[1]: ', Tlist[1])

    nlist = [] #number of columns of each matrix in Rlist (it may not always be the same value for all matrices)
    for k in Rlist:
        nlist.append(len(k.mat[0]))
    
    p = np.zeros(N).tolist()
    for k in range(shots):
        # generate uniformly random bitstrings u
        ulist = []
        for l in range(len(Rlist)):
            tmp = np.random.choice([0,1],nlist[l])
            print('tmp = ', tmp)
            ulist.append(tmp.tolist())
        #sample the probability distribution with the generated bitstrings
        samp = sample_shot(Rlist, Tlist, clist, ulist)
        #add to the sampled distribution
        p = add_vec(p, samp)
        print('sample = ' ,samp)
        print('p = ',p)
    return p

#============================================================

#associate coeffs to keys, according to the first element of dlist
# c = {}
# for k in coeffs:
#     c[key[k]] = k

#============================================================

#coefficient in the linear combination
# coeff = [4/np.sqrt(17), np.pi*1j/np.sqrt(34)] #first element corresponds to stab[0], etc
# nrm = 2 + pow(np.pi, 2)/16 - pow(np.pi, 4)/256
# nrm = 2 + pow(np.pi, 2)/16 - pow(np.pi, 4)/256 + pow(np.pi, 6)/4096
# coeff = [np.sqrt(2)/nrm, (np.pi*1j/4- pow(np.pi, 2)/16)/(nrm)] #TODO: the problem is probably also here, but one thing is sure it's that there must be not abs here 

#============================================================

#This codes comes from the function "project" in weaksim.py
#deleting the columns l where R[qb][l] = 1 
        # for k in sorted(idx, reverse=True):
        #     if len(R.mat[0]) > 1: #more than one column
        #         R.del_col(k)
        # if len(idx) == 1:
        #     idx = idx[0]
        #     #TODO: faire comme pour le else
        # else: 
        #     # print('Big nope')
        #     # idx = idx[0]
        # if len(R.mat[0]) > 1: #more than one column
        #     print('bim')
        #     print(len(R.mat[0]))
        #     R.del_col(idx) 

#============================================================

def sample_shot(Rlist, Tlist, coefflist, ulist): #TODO: corriger, le problème vient du val += ... probablement -> normalement c'est corrigé
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
        plist.mat[k] = sample_shot(Rlist, Tlist, coefflist, ulist)
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
    pfin = add_dict(drac, coefflist, True)
    return pfin 

#============================================================

def flatten(mat):
    #Flattens the matrix row wise, so that the end matrix has only 1 row
    out = []
    for k in range(len(mat)):
        for l in mat[k]:
            out.append(l)
    return out

#============================================================
