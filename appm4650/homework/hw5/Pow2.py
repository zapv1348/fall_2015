import numpy as np

A=np.array([[2,-1,1],[-1,3,2,],[1,2,3]])
x0=np.array([[1],[1],[1]])
tol=10**(-4)


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


