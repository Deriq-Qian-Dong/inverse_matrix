#coding=utf-8
import math
def mul(A,B):
    N=int(math.sqrt(len(A)))
    C=[0]*len(A)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i*N+j]+= A[i*N + k] * B[k*N + j]
    for i in range(len(A)):
        if(abs(C[i])<10**-10):
            C[i]=0
    return  C

def LUP_Descomposition(A,L,U,P):
    row = 0
    N = int(math.sqrt(len(A)))
    for i in range(N):
        P[i]=i
    for i in range(N-1):
        p=0.0
        for j in range(i,N):
            if(abs(A[j*N + i])>p):
                p = abs(A[j * N + i])
                row = j
        if(p==0):
            print("矩阵奇异，无法计算逆")
            for a in A:
                print a
            return
        tmp=P[i]
        P[i]=P[row]
        P[row]=tmp
        tmp2=0.0
        for j in range(N):
            tmp2=A[i*N + j]
            A[i * N + j] = A[row * N + j]
            A[row * N + j] = tmp2
        u = A[i * N + i]
        l = 0.0
        for j in range(i+1,N):
            l = A[j * N + i] / u
            A[j * N + i] = l
            for k in range(i+1,N):
                A[j * N + k] = A[j * N + k] - A[i * N + k] * l
    for i in range(N):
        for j in range(i):
            if(i!=j):
                L[i * N + j] = A[i * N + j]
            else:
                L[i * N + j] = 1
        for k in range(i,N):
            U[i * N + k] = A[i * N + k]
    return A,L,U,P
def LUP_Solve(L,U,P,b):
    N = int(math.sqrt(len(L)))
    x=[0.0]*N
    y=[0.0]*N
    for i in range(N):
        y[i] = b[P[i]]
        for j in range(i):
            y[i] = y[i] - L[i * N + j] * y[j]

    for i in range(N-1,-1,-1):
        x[i]=y[i]
        for j in range(N-1,i,-1):
            x[i] = x[i] - U[i * N + j] * x[j]
        x[i] /= U[i * N + i]
    return x

def getNext(i,m,n):
    return (i%n)*m + i / n
def getPre(i,m,n):
    return (i % m) * n + i / m
def movedata(mtx,i,m,n):
    temp=mtx[i]
    cur=i
    pre=getPre(cur,m,n)
    while(pre!=i):
        mtx[cur] = mtx[pre]
        cur = pre
        pre = getPre(cur, m, n)
    mtx[cur]=temp
    return mtx

def transpose(mtx,m,n):
    for i in range(m*n):
        next=getNext(i,m,n)
        while(next>i):
            next = getNext(next, m, n)
        if(next==i):
            mtx=movedata(mtx, i, m, n)
    return mtx

def LUP_solve_inverse(A):
    A=[float(a) for a in A]
    N=int(math.sqrt(len(A)))
    A_mirror=[0.0]*len(A)
    inv_A=[0.0]*len(A)
    inv_A_each=[0.0]*N
    b=[0.0]*N
    for i in range(N):
        L=[0.0]*len(A)
        U=[0.0]*len(A)
        P=[0]*N
        for j in range(N):
            b[j]=0
        b[i]=1
        A_mirror=[a for a in A]
        A_mirror, L, U, P=LUP_Descomposition(A_mirror, L, U, P)
        inv_A_each = LUP_Solve(L, U, P, b)
        inv_A[i*N:i*N+N]=[a for a in inv_A_each]
    inv_A=transpose(inv_A, N, N)
    return  inv_A
if __name__=="__main__":
    A=[1,2,3,1,0,-1,0,1,1]
    invOfA = LUP_solve_inverse(A)
    for a in invOfA:
        print a

