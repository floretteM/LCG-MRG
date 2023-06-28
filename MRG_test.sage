# Please load MRG.sage

import time

#################
### a,N known ###
#################

### CVP ###

def test_CVP(n,k,ell,m):
    res = 0
    for round in range(100):
        N = ZZ.random_element(2^(n-1),2^n-1)
        a = []
        for i in range(k):
            a.append(ZZ.random_element(1,N-1))
        x_init,Y,x = MRG(k,a,N,m,ell,[])
        H = transform(Y,ell)
        x_init_bis = attack_CVP(k,H,a,N,m,ell)
        if x_init_bis == x_init:
            res+=1
    return(RR(res/100))

def test_time_CVP(n,k,ell,m):
    res = 0
    instances = []
    for round in range(100):
        N = ZZ.random_element(2^(n-1),2^n-1)
        a = []
        for i in range(k):
            a.append(ZZ.random_element(1,N-1))
        x_init,Y,x = MRG(k,a,N,m,ell,[])
        H = transform(Y,ell)
        instances.append([N,a,x_init,H])
    start = time.time()
    for round in range(100):  
        inst = instances[round]
        N = inst[0]
        a = inst[1]
        x_init = inst[2]
        H = inst[3]
        x3 = attack_CVP(k,H,a,N,m,ell)
    stop = time.time()
    return(RR((stop-start)/100))

### Frieze ###
def test_frieze(n,k,ell,m):
    res = 0
    for round in range(100):
        N = ZZ.random_element(2^(n-1),2^n-1)
        a = []
        for i in range(k):
            a.append(ZZ.random_element(1,N-1))
        x_init,Y,x = MRG(k,a,N,m,ell,[])
        H = transform(Y,ell)
        x_init_bis = attack_frieze(k,H,a,N,m,ell)
        if x_init_bis == x_init:
            res+=1
    return(RR(res/100))

def test_time_frieze(n,k,ell,m):
    res = 0
    instances = []
    for round in range(100):
        N = ZZ.random_element(2^(n-1),2^n-1)
        a = []
        for i in range(k):
            a.append(ZZ.random_element(1,N-1))
        x_init,Y,x = MRG(k,a,N,m,ell,[])
        H = transform(Y,ell)
        instances.append([N,a,x_init,H])
    start = time.time()
    for round in range(100):  
        inst = instances[round]
        N = inst[0]
        a = inst[1]
        x_init = inst[2]
        H = inst[3]
        x3 = attack_frieze(k,H,a,N,m,ell)
    stop = time.time()
    return(RR((stop-start)/100))


#########################
### a unkown, N known ###
#########################

### Stern simple ###

def find_prime(n):
    N=ZZ.random_element(2^(n-1),2^n-1)
    while N.is_prime == False:
        N=ZZ.random_element(2^(n-1),2^n-1) 
    return(N)

def test_stern_simple(n,k,ell,m,r):
    res = 0
    for round in range(100):
        N = find_prime(n)
        a = []
        for i in range(k):
            a.append(ZZ.random_element(1,N-1))
        x_init,Y,x = MRG(k,a,N,m,ell,[])
        H = transform(Y,ell)
        a_0 = attack_stern_simple(k,H,N,r,ell)
        if a_0 == a[0]:
            res+=1
    return(RR(res/100))

def test_time_stern_simple(n,k,ell,m,r):
    res = 0
    instances = []
    for round in range(100):
        N = find_prime(n)
        a = []
        for i in range(k):
            a.append(ZZ.random_element(1,N-1))
        x_init,Y,x = MRG(k,a,N,m,ell,[])
        H = transform(Y,ell)
        instances.append([N,a,x_init,H])
    start = time.time()
    for round in range(100):  
        inst = instances[round]
        N = inst[0]
        a = inst[1]
        x_init = inst[2]
        H = inst[3]
        a_0 = attack_stern_simple(k,H,N,r,ell)
    stop = time.time()
    return(RR((stop-start)/100))

def test_full_stern(n,k,a,N,m,ell,d,r):
    x_init = []
    x_init,Y,X = MRG(k,a,N,m,ell)
    coeff = find_coeff_full(n,Y,d,r,ell)
    tmp = sum([coeff[i]*X[i] for i in range(d)])
    print(':p', tmp %N)
    list_P = find_poly_full(n,k,Y,d,r,ell)
    for P in list_P:
        print(P,P(a) % N)