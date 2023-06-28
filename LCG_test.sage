# Please load LCG.sage
from time import *

#################
### a,N known ###
#################

### CVP ###

def test_CVP(n,ell,m):
    res = 0
    instances = []
    for round in range(100):
        N = ZZ.random_element(2^(n-1),2^n-1)
        a = ZZ.random_element(1,N-1)
        x0=0
        x0,Y = LCG(a,N,m,ell,x0)
        H = transform(Y,ell)
        instances.append([N,a,x0,H])
    T = time()
    for i in range(100):  
        [N,a,x0,H] = instances[i]
        x0bis = attack_CVP(H,a,N,m,ell)
        if x0bis == x:
            res+=1
    T = time()- T
    return(1.0*res/100,T/100)

### Frieze ###

def test_frieze(n,ell,m):
    res = 0
    instances = []
    for round in range(100):
        N = ZZ.random_element(2^(n-1),2^n-1)
        a = ZZ.random_element(1,N-1)
        x0=0
        x0,Y = LCG(a,N,m,ell,x0)
        H = transform(Y,ell)
        instances.append([N,a,x0,H])
    T = time()
    for i in range(100):  
        [N,a,x0,H] = instances[i]
        x0bis = attack_frieze(H,a,N,m,ell)
        if x0bis == x:
            res+=1
    T = time()- T
    return(1.0*res/100,T/100)

#########################
### a unkown, N known ###
#########################

### Stern simple ###

def find_prime(n):
    N=ZZ.random_element(2^(n-1),2^n-1)
    while N.is_prime() == False:
        N=ZZ.random_element(2^(n-1),2^n-1) 
    return(N)

def generate_prime(n):
    list_N=[]
    for i in range(100):
        list_N.append(find_prime(n))
    return(list_N)

def test_polynomial_simple(n,ell,mp):
    d = (mp+1)//2
    instances=[]
    for i in range(100):
        N = ZZ.random_element(2^(n-1),2^n-1)
        a = ZZ.random_element(1,N-1)
        x0=0
        x0,Y = LCG(a,N,mp,ell,x0)
        H = transform(Y,ell)
        instances.append([x0,a,H,N])
    print('end of initialization phase')
    res=0
    T = time()
    for i in range(100):
        [x0,a,H,N] = instances.pop()
        coeff = find_polynomial_simple(H,N,ell)
        P = sum([coeff[i]*a^i for i in range(d)])
        if P % N ==0:
            res+=1
    T=time()-T
    return(1.0*res/100,T/100)    

def test_simple_stern(n,ell,mp,list_N):
    d=(mp+1)//2
    m=d+mp
    instances=[]
    for i in range(100):
        N = list_N[i]
        a = ZZ.random_element(1,N-1)
        x0=0
        x0,Y = LCG(a,N,m,ell,x0)
        H = transform(Y,ell)
        instances.append([x0,a,H,N])
    print('end of initialization phase')
    res=0
    T = time()
    for i in range(100):
        [x0,a,H,N] = instances.pop()
        temp = attack_simple_stern(H,N,mp,ell)
        if temp != None:
            [x0bis,abis] = temp
            if ZZ(abis) == a and x0bis == x0:
                res+=1
    T=time()-T
    return(1.0*res/100,T/100)

def test_simple_stern_alt(n,ell,m):
    d = (m+1)//2
    N=2^n
    instances=[]
    for i in range(100):
        a = ZZ.random_element(1,N-1)
        x0,Y = LCG(a,N,m,ell,0)
        H = transform(Y,ell)
        instances.append([x0,a,H,N])
    print('end of initialization phase')
    res=0
    T = time()
    for i in range(100):
        [x0,a,H,N] = instances.pop()
        temp = attack_simple_stern_alt(H,n,ell)
        if temp != None:
            x0bis,abis = temp
            if ZZ(abis) == a and x0bis == x0:
                res+=1
    T=time()-T
    return(1.0*res/100,T/100)

def find_theo_ell(n,mp):
    d = (mp+1)//2
    r = (mp+1) - d
    ell = 1
    a = sqrt(d+r*d^2)*2^(ell-1)*2^(n/d)
    b = sqrt(r+d)*(2^(d*(ell-1))*2^(n*r))^(1/(d+r))
    while a <= b :
        ell +=1
        a = sqrt(d+r*d^2)*2^(ell-1)*2^(n/d)
        b = sqrt(r+d)*(2^(d*(ell-1))*2^(n*r))^(1/(d+r))
    return(ell-1)
### Coppersmith ###

def test_coppersmith(n,ell,m):
    res = 0
    instances = []
    for round in range(100):
        N = ZZ.random_element(2^(n-1),2^n-1)
        a = ZZ.random_element(1,N-1)
        x0=0
        x0,Y = LCG(a,N,m,ell,x0)
        H = transform(Y,ell)
        instances.append([N,a,x0,H])
    start = time()
    for i in range(100):  
        [x0,a,H,N] = instances.pop()
        x0bis = attack_Coppersmith(H,N,ell)
        if x0bis == x0:
            res+=1
    stop = time()
    return(1.0*res/100,(stop-start)/100)

###########################
### a unkown, N unknown ###
###########################

### Stern ###

def test_find_coeff_full(n,ell,r,list_N):
    res = 0
    if n> 100:
        d = ceil(sqrt(2*(n-ell)*r))
    else :
        d = ceil(sqrt(2*n*r))
    m = r+d
    instances = []
    for i in range(100):
        N = list_N[i]
        a = ZZ.random_element(1,N-1)
        x0=0
        x0,Y = LCG(a,N,m,ell,x0)
        instances.append([N,a,x0,Y])
    T=time()
    for i in range(100):
        [N,a,x0,Y] = instances.pop()
        b=(n-ell +log(d,2) -1)*r/(d-r)
        B=2^b
        K = ceil(sqrt(d)*2^(d-1)/2*B)
        temp = find_coeff_full(Y,K,r,d,ell)
        P = sum([temp[i]*a^i for i in range(d)])
        if P%N == 0:
            res+=1
    T=time()-T
    return(T/100,RR(res/100))

def test_stern(n,ell,r,list_N):
    res = 0
    if n> 100:
        d = ceil(sqrt(2*(n-ell)*r))
    else :
        d = ceil(sqrt(2*n*r))
    m = 3*r+d
    instances = []
    for i in range(100):
        N = list_N[i]
        a = ZZ.random_element(1,N-1)
        x0=0
        x0,Y = LCG(a,N,m,ell,x0)
        instances.append([N,a,x0,Y])
    T=time()
    for i in range(100):
        [N,a,x0,Y] = instances.pop()
        temp = attack_stern(n,Y,r,d,ell)
        if temp != None:
            x0bis,abis,Nbis = temp
            if x0bis == x0 and ZZ(abis) == a and Nbis == N:
                res+=1
    T=time()-T
    return(T/100,RR(res/100))
