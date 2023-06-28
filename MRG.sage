from fpylll import *
import time

""" a: list_multiplier
    N: modulus
    m: nb of outputs
    ell: nb of truncated bits
"""

def MRG(k,a,N,m,ell,x_init = []):
    if x_init == []:
        for i in range(k):
            x_init.append(ZZ.random_element(1,N-1))
    x = [int(xj) for xj in x_init]
    for j in range(k,m):
        xj = sum([a[i]*x[j-k+i] for i in range(k)]) % N
        x.append(xj)
    res=[(xj>>ell) for xj in x ]
    return(x_init,res,x)

def transform(Y,ell):
    if ell == 0:
        return Y
    H=[]
    for y in Y:
        H.append(y*2^ell + 2^(ell-1))
    return(H)

#################
### a,N known ###
#################

### CVP ###

def attack_CVP(k,H,a,N,m,ell):
    #straight forward attack using a CVP solver. Inconvenient: expensive

    L = matrix(m)
    for j in range(m):
        if j<k:
            L[j,j] = 1
        else:
            L[j,j] = N
            for i in range(k):
                L[i,j] = sum([a[s]*L[i,j-k+s] for s in range(max(0,k-j),k)])
    L  = IntegerMatrix.from_matrix(L)
    L = fplll.lll.lll_reduction(L)
    target = tuple(vector(H[:m]))
    x_init = CVP.closest_vector(L,target)[:k]

    return list(x_init)

### Frieze ###

def attack_frieze(k,H,a,N,m,ell):

    A = matrix(m)
    for j in range(m):
        if j<k:
            A[j,j] = 1
        else:
            A[j,j] = -1
            for i in range(k):
                A[i,j] = sum([a[s]*A[i,j-k+s] for s in range(max(0,k-j),k)]) % N
    for j in range(k):
        A[j,j] = N
    A= A.transpose()
    #print('A = ',A)
    A = IntegerMatrix.from_matrix(A)
    A = fplll.lll.lll_reduction(A)
    #print('A.LLL=',A)
    A = matrix(A)
    v = vector(H[:m])
    temp = (A*v)/N
    K = vector(ZZ,[round(t) for t in temp])
    X = A.solve_right(N*K)
    x_init = X[:k]
    return list(x_init)

#########################
### a unkown, N known ###
#########################    

### Stern simple ###

def my_gcd(p1,p2,N): 
    #for some reasons, the gcd does not work anymore in this case with the new version of sage so I had to made my own
    R=p1.parent()
    a0 = var('a0')
    a0=R(a0)
    r1 = p1
    r2 = p2
    if r2.degree()>r1.degree():
        r1,r2 = r2,r1
    while r2 != 0:
        r3 = r1
        q=0
        while r3.degree() >= r2.degree():
            d = r3.degree() - r2.degree()
            a = r3.leading_coefficient()
            b = r2.leading_coefficient()
            if b.is_unit() == False:
                return(r3) #we shortcut the GCD in case of a bad situation
            else :
                r3 = r3 - a*b^(-1)*a0^d*r2
                q += a*b^(-1)*a0^d
        assert(r1 == q*r2+r3)
        r1 = r2
        r2 = r3
    return(r1)

def find_coeff_simple(H,N,r,ell):
    """ Input:  H is the modified ouputs of a MRG (public)
                N he modulus (public)
                ell the number of truncated bits  (public)
                r is a parameter to choose

                
        Output: a vector "target" such that <target,H_i> = 0 [N]. """

    block_H = matrix(r)
    for i in range(r):
        for j in range(r):
            block_H[i,j] = H[i+j]
    M=block_matrix(2,2,[2**(ell-1)*identity_matrix(r),block_H,matrix(r),N*identity_matrix(r)])
    M = IntegerMatrix.from_matrix(M)
    M=fplll.lll.lll_reduction(M)
    M=matrix(M)
    target = vector(ZZ,M[0][:r]/(2**(ell-1)))
    return(target)

def cleaning(Q):
    i=0
    while Q[i] == 0:
        i+=1
    return(i)
    
def find_poly_simple(k,H,N,r,ell):
    """ Input:  H is the modified ouputs of a MRG (public)
                N he modulus (public)
                ell the number of truncated bits  (public)
                k is the order of the recursive sequence
                r is a parameter to choose

                
        Output: Q a polynomial such that Q(a_0) = 0 [N] """

    coefs = find_coeff_simple(H,N,r,ell)

    list_var = var(['a'+ str(i) for i in range(k)])
    R = PolynomialRing(ZZ,list_var)
    a = [R(variable) for variable in list_var]

    L = matrix(R,k,r) #this matrix will contain the monomials
    for j in range(r):
        if j<k:
            L[j,j] = 1
        else:
            for i in range(k):
                L[i,j] = sum([a[s]*L[i,j-k+s] for s in range(max(0,k-j),k)])  
  
    list_P = []
    for i in range(k):
        Pi = sum([coefs[j]*L[i,j] for j in range(r)])
        list_P.append(Pi)

    for i in range(k-1,0,-1):
        list_temp = [] 
        P0 = list_P[0]  
        for j in range(1,i+1):
            Pj = list_P[j]
            list_temp.append(P0.resultant(Pj,a[i]))
        list_P = list_temp
    Q=list_P[0]

    R2.<a0> = Integers(N)[]
    Q=R2(Q)
    """
    i = cleaning(Q)
    Q=Q/a0^i
    Q=R2(Q)
    """
    return(R2(Q))


def attack_stern_simple(k,H,N,r,ell):
    """ Input: H is the modified ouputs of a LCG (public)
                N he modulus (public)
                ell the number of truncated bits  (public)
                k is the order of the recursive sequence
                r is a parameter to choose

                
        Output: a_0 """ 
    #only with N prime otherwise resultant not implemented
    i=0 
    #R.<a0> = Integers(N)[]
    Qi = find_poly_simple(k,H,N,r,ell)
    Q=Qi
    #return(Q)
    while Q.degree() > 1 and i<10:
        i+=1
        Qi = find_poly_simple(k,H[i:],N,r,ell)
        
        Q = my_gcd(Q,Qi,N)

    if Q.degree() > 1:
        return(0)

    if Q[1] == 0:
        return 0
    a_0 = -Q[0]/Q[1] 
    return(a_0)


#######################
### a and N unknown ###
#######################

### NOT CONVINCING ###

def find_coeff_full(n,Y,d,r,ell):
    block_Y = matrix(d,r)
    b=(n-ell +log(d,2) -1)*r/(d-r)
    B=2^b
    K = ceil(sqrt(d)*2^(d-1)/2*B)
    for i in range(d):
        for j in range(r):
            block_Y[i,j] = Y[i+j]
    M=block_matrix(1,2,[identity_matrix(d),K*block_Y])
    M = IntegerMatrix.from_matrix(M)
    M=fplll.lll.lll_reduction(M)
    M=matrix(M)
    if M[0][d:] == 0:
        target = vector(ZZ,M[0][:d])
        return(target)
    else :
        print("sum(mu_i H_i) != 0")

def find_poly_full(n,k,Y,d,r,ell):
    coefs = find_coeff_full(n,Y,d,r,ell)
    list_var = var(['a'+ str(i) for i in range(k)])
    R = PolynomialRing(ZZ,list_var)
    a = [R(variable) for variable in list_var]

    L = matrix(R,k,r) #this matrix will contain the monomials
    for j in range(r):
        if j<k:
            L[j,j] = 1
        else:
            for i in range(k):
                L[i,j] = sum([a[s]*L[i,j-k+s] for s in range(max(0,k-j),k)])  
  
    list_P = []
    for i in range(k):
        Pi = sum([coefs[j]*L[i,j] for j in range(r)])
        list_P.append(Pi)

    return(list_P)

    for i in range(k-1,0,-1):
        list_temp = [] 
        P0 = list_P[0]  
        for j in range(1,i+1):
            Pj = list_P[j]
            list_temp.append(P0.resultant(Pj,a[i]))
        list_P = list_temp
    Q=list_P[0]

    R2.<a0> = ZZ[]
    Q=R2(Q)
    return(Q)


    """
    i = cleaning(Q)
    Q=Q/a0^i
    Q=R2(Q)
    return(R2(Q))"""


def attack_stern_full(n,k,H,d,r,ell):
    #d = ceil(sqrt(2*(n-ell)*r))
    i=0
    #R.<a0> = Integers(N)[]
    Qi = find_poly_full(n,k,H,d,r,ell)
    Q=Qi
    while Q.degree() > 1 and i<10:
        i+=1
        Qi = find_poly_full(k,H[i:],d,r,ell)
        print("Q",Q)
        print("Qi",Qi)
        Q = Q.resultant(Qi,a0)
        
        j =cleaning(Q)
        if j!=0:
            print('cleaning needed')
    
    if Q.degree() > 1:
        return(0)
    #print('i',i)
    if Q[1] == 0:
        return 0
    a_0 = -Q[0]/Q[1] 
    return(a_0)

    return(Q)










