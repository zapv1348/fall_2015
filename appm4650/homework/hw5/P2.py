import numpy as np


A=np.array([[3,0,1],[0,5,0],[-1,1,-1]])
n=A[:,1].size
x1=np.zeros_like(A)
x0=np.array([[0.5,-0.1,0.4],[0,0.2,0],[-0.4,0.3,-1.5]])
I=np.identity(3)

x1=np.dot(x0,(2*I-np.dot(A,x0)))

print("x1 is")
print(x1)

Ain=np.linalg.inv(A)

print("A^(-1)-x0")
print(Ain-x0)

print("A^(-1)-x1")
print(Ain-x1)
