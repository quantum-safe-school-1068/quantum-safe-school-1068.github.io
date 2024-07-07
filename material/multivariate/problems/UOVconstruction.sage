# Multivariate Cryptography
# Introduction to Quantum-Safe Cryptography (IBM Zurich)
# Simona Samardjiska
# Programming assignment


# Implementation of the UOV signature scheme
# the implementation matches the general design of UOV
# details not relevant to the design paradigm are neglected

######################################
# helper

# evaluate Multivariate map 
def MQeval(F,x,y):
    return [(x.transpose()*M*y)[0,0] for M in F]

def MQevalleft(F,x):
    return [ x.transpose()  *M for M in F]

#######################################
# UOV Key generation

# generate random invertible matrix 
def RandomInvertible(n):
    retM = Matrix(K,n,n)
    while not retM.is_invertible():
        retM = Matrix(K,[[K.random_element() for j in range(n)] for i in range(n)])
    return retM

# turn random matrix to upper diagonal
def SquareToUpper(M):
    n = M.ncols()
    retM = M
    for i in range(n):
        for j in range(i+1,n):
            retM[i,j] += retM[j,i]
            retM[j,i] = K(0) 
    return retM
            
# turn upper diagonal matrix to symmetric
def UpperToSymmetric(M):
    return M+M.transpose() 

# generates a UOV public-private key pair
# we neglect the 0s in the public/private key matrices when calculating key sizes, 
# and assume a triangular form is equivalent to a compressed representation 
def Keygen(q,n,m):
    Central_Map = [];
    for k in range(m):
        P  = Matrix(K,n,n)
        for i in range(n-m):
            for j in range(i,n):
                P[i,j] = K.random_element()
        Central_Map.append(P)
    S = RandomInvertible(n)
    Public_Key = [SquareToUpper(S.transpose()*M*S) for M in Central_Map]
    
    return Public_Key, S, Central_Map

def Sign(hash, S, Central_Map):
    Central_vinegar = [M.submatrix(0,0,n-m,n-m) for M in Central_Map]
    Central_oilvinegar = [M.submatrix(0,n-m,n-m,m) for M in Central_Map]
    coef_matrix = Matrix(K,m,m)
    while not coef_matrix.is_invertible():
        vinegar_vector=Matrix(K,n-m,1,[K.random_element() for j in range(n-m)])
        coef_matrix = Matrix(K,[(vinegar_vector.transpose() * Central_oilvinegar[i]).list()  for i in range(m)])
    coef_vector = Matrix(K,m,1,[(hash[i,0] - vinegar_vector.transpose() * Central_vinegar[i] * vinegar_vector)[0,0]  for i in range(m)])
    oil_vector = coef_matrix.inverse() * coef_vector
    signature = S.inverse() * vinegar_vector.stack(oil_vector)
    return signature

def Verify(hash, signature, Public_Key):
    hashprime = MQeval(Public_Key,signature,signature)
    return hash.list() == hashprime


