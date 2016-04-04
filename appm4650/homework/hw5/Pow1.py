import numpy as np

A=np.array([[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]])
x0=np.array([[1],[1],[1],[1]])
x1=np.array([[1],[1],[5],[1]])
tol=10**(-6)


xold=np.zeros_like(x0)
n=0
while (np.abs(xold-x0)>tol).all():
    xold=x0
    s=np.dot(A,x0)
    x0=s/np.linalg.norm(s)
    xt=np.transpose(x0)
    lam=np.dot(xt,s)/(np.dot(xt,x0))
    s="iteration=%d\n"% n
    print(s)
    print(x0)
    print(lam)
    n=n+1


print()
n=0
while (np.abs(xold-x1)>tol).all():
    xold=x1
    s=np.dot(A,x1)
    x1=s/np.linalg.norm(s)
    xt=np.transpose(x1)
    lam=np.dot(xt,s)/(np.dot(xt,x1))
    s="iteration=%d\n"% n
    print(s)
    print(x1)
    print(lam)
    n=n+1
