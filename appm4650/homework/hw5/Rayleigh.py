import numpy as np

A=np.array([[2,-1,1],[-1,3,2,],[1,2,3]])
x0=np.array([[1],[1],[1]])
lam=10
I=np.identity(3)

n=0
while n<4:
    f=np.linalg.inv(A-lam*I)
    s=np.dot(f,x0)
    xt=np.transpose(x0)
    lam=np.dot(xt,np.dot(A,x0))/np.dot(xt,x0)
    x0=s/np.linalg.norm(s)
    r=("iteration= %d\n")%(n)
    print(r)
    print("x is:")
    print(x0)
    print("the eigenvalue is:")
    print(lam)
    n=n+1


