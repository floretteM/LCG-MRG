from fpylll import *
import time
###################
### Coppersmith ###
###################

def coppersmith(G,monomials,N,l,m):
    """ retrieve a little common root of polynomials over a certain modulo
    
    Input:  G a list of polynomials
                monomials the list of monomials appearing in G (in a certain order)
                N the moduli
                l the size of the unknown root
                m the multiplicity of the root in the polynomials of G"""

    nb_poly = len(G)
    nb_mono = len(monomials)

    # Construction of the matrix described in section 2
    
    M = matrix(QQ,nb_mono+nb_poly)
    for i in range(nb_mono):
        if monomials[i] == 1: 
            M[i,i] = 1
        else:
            M[i,i] = (2^l)^(-1*monomials[i].degree())

    for j in range(nb_poly):
        M[nb_mono+j,nb_mono+j] = N^m
        for i in range(nb_mono):
            M[i,nb_mono+j] = G[j].monomial_coefficient(monomials[i])
                   
    # LLL
    M=M.LLL()

    
    
    # recuperation of the smallest vector with the correct form
    
    i=0
    
    while i < nb_mono+nb_poly and  abs(M[i,0]) == 0:
        i=i+1
    if i != nb_mono+nb_poly : 
        v = M[i,0]^-1*M[i]
        v_final = 2^l*v

        return v_final
        
    return 0 #In the case the algorithm fails



m1 = 2^32 - 209
m2 = 2^32-22853

a11 = 0
a12 = 1403580
a13 = 810728

a21 = 527612
a22 = 0
a23 = 1370589    
    
    
def MGR(a11,a12,a13,a21,a22,a23,m1,m2,nb_out):

    A = CRT(a11,a21,m1,m2)
    B = CRT(a12,a22,m1,m2)
    C = CRT(a13,a23,m1,m2)
    
    X0 = ZZ.random_element(0,m1*m2)
    X1 = ZZ.random_element(0,m1*m2)
    X2 = ZZ.random_element(0,m1*m2)
    
    X = [X0,X1,X2]
    
    while len(X) < nb_out:
        new_el = (A*X[-1] + B*X[-2] + C*X[-3]) %(m1*m2)
        X.append(new_el)
        
    x = [Xi % m1 for Xi in X]
    y = [Xi % m2 for Xi in X]
    
    z=[]
    for i in range(nb_out):
        zi = (x[i]-y[i]) % m1
        z.append(zi)
        
    return X,x,y,z
    
def attack(a11,a12,a13,a21,a22,a23,m1,m2,n,z):
    """ return the seed X0, X1, X2 of the MRG32
    If the attack fails, return nothing."""

    #start = time.time()
    
    u,v = xgcd(m1,m2)[1:3]
    A = CRT(a11,a21,m1,m2)
    B = CRT(a12,a22,m1,m2)
    C = CRT(a13,a23,m1,m2)    
    
    k = []

    m=min(len(z)-1,7) 
    #We do not need more than 8 outputs to attack so we put a limit there.
    for r in range(2^m):
        sign = bin(r)[2:]
        sign = '0'*(m-len(sign))+sign
        for i in range(m):
            ziprime = z[i] - int(sign[i])*m1
            ki = -u*ziprime % m2
            k.append(ki)
        R.<v0,v1,v2,v3,v4,v5,v6> = ZZ[]
        v = [v0,v1,v2,v3,v4,v5,v6]
        G=[]
        for i in range(3,m):
            G.append( k[i]*m1+v[i] - A*(k[i-1]*m1 + v[i-1]) - B*(k[i-2]*m1+v[i-2]) - C*(k[i-3]*m1+v[i-3]))
        monomes = [v0^0,v0,v1,v2,v3,v4,v5,v6]
        res = coppersmith(G,monomes,m1*m2,n,1)[1:m+1]
        X_0 = k[0]*m1+res[0]
        X_1 = k[1]*m1+res[1]
        X_2 = k[2]*m1+res[2]
        X=[X_0,X_1,X_2]
        while len(X) < m+1:
            new_el = (A*X[-1] + B*X[-2] + C*X[-3]) %(m1*m2)
            X.append(new_el)
        x_f = X[-1] %m1
        y_f = X[-1] %m2
        z_f = (x_f-y_f) %m1
        if z_f == z[-1] :  
            #stop = time.time()-start
            #print(stop)   
            return(X)
    
def test_MGR(n,nb_out):
    m1 = 0
    while m1.is_prime() == False:
        m1 = ZZ.random_element(0,2^n-1)
        
    m2 = 0
    while m2.is_prime() == False:
        m2 = ZZ.random_element(0,2^n-1)        
        
    a11 = ZZ.random_element(0,2^n-1)
    a12 = ZZ.random_element(0,2^n-1)
    a13 = ZZ.random_element(0,2^n-1)

    a21 = ZZ.random_element(0,2^n-1)
    a22 = ZZ.random_element(0,2^n-1)
    a23 = ZZ.random_element(0,2^n-1)
    
    X,x,y,z = MGR(a11,a12,a13,a21,a22,a23,m1,m2,nb_out)
    Xbis = attack(a11,a12,a13,a21,a22,a23,m1,m2,n,z)
    if X==Xbis:
        return True
    else:
        return False
