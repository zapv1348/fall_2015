import numpy as np

#dat cofactor method for inverse matrix finding
#inverse is 1/det*adjugate which is checkerboard of +,- on the og matrix transposed.

A=np.array([[3,0,1],[0,5,0],[-1,1,-1]])
n=A[:,1].size
Ain=np.zeros_like(A)


def twobytwo(A,i,j):
    a=(i+1)%n
    b=(i+2)%n
    c=(j+1)%n
    d=(j+2)%n
    if (i+j)%2==0:
        return A[a,c]*A[b,d]-A[a,d]*A[b,c]
    else:
        return A[a,d]*A[b,c]-A[a,c]*A[b,d]

B=np.zeros_like(A)
for i in range(n):
    for j in range(n):
        B[i,j]=twobytwo(A,i,j)

#find the adjunct matrix
for i in range(n):
    for j in range(n):
        if (i+j)%2==0:
            Ain[j,i]=B[i,j]
        else:
            Ain[j,i]=-B[i,j]

det=B[0,0]*A[0,0]-B[0,1]*A[0,1]+B[0,2]*A[0,2]
Ain=1/det*Ain

print("A inverse is:")
print(Ain)

