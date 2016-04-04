import numpy as np

h=0.25
a=0
b=2
x=0;
est=0;


def my_range(start,end,step):
    while start<=end:
        yield start
        start+=step

def exp_thing(z):
    return np.exp(-(x**2))*x**2

for x in my_range(a+h,b-h,2*h):
    est=h/3*(exp_thing(x-h)+4*exp_thing(x)+exp_thing(x+h))+est

print("estimate")
print(est)
