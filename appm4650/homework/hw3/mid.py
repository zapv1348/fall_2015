import numpy as np

h=0.25
a=0
b=2
x1=0
x=0
est=0;


def my_range(start,end,step):
    while start<=end:
        yield start
        start+=step


for x in my_range(a+h/2,b-h/2,h):
    est=np.exp(-(x**2))*x**2+est

s=a+h/2
r=b-h/2
n=(int)((r-s)/h)
est=(b-a)*est/n
print("delta x:")
print((b-a)/n)
print("estimate")
print(est)
