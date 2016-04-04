import numpy as np

h=0.25
a=0
b=2
x=0;
x1=a;
est=0;


def my_range(start,end,step):
    while start<=end:
        yield start
        start+=step


for x in my_range(a+h,b,h):
    est=np.exp(-(x**2))*x**2+np.exp(-(x1**2))*x1**2+est
    x1=x

est=h*est/2
print("estimate")
print(est)
