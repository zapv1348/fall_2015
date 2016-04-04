import numpy as np

A=np.array([[4., 3., 0.],[3., 4., -1.],[0., -1., 4.]])
b=np.array([24.,30.,-24.])

n=A.size

xold=np.array([2.,2.,2.])

x=np.zeros_like(xold);
j=0
w=1.25


while ~(np.allclose(xold,x,rtol=1e-14)):
    xold=x
    print ("current val:",xold)
    x=np.zeros_like(xold)

    for i in range(A.shape[0]):
        a1=np.dot(A[i,:i],x[:i])
        a2=np.dot(A[i,i+1:],xold[i+1:])
        s=(b[i]-a1-a2)/A[i,i]
        x[i]=(1-w)*xold[i]+w*s
    print(j)
    j=j+1

print("Solution")
print(xold)
print("b:")
print(np.dot(A,xold))
