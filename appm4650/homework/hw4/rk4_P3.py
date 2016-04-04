import numpy as np

ydiff="x+y"

y0=0.0
x0=0.0

xmax=0.5
h=0.1

ya=np.array([y0, 0, 0, 0, 0, 0])
xa=np.array([x0,h,2*h,3*h,4*h,5*h])


def my_range(start,end,step):
    while start<=end:
        yield start
        start+=step

r="iteration |       x      |       y        |"
adswnloadswnloads
print(r)


for x in my_range(x0,xmax-h,h):
    n=x*10
    y=ya[n]
    x=xa[n]
    k1=eval(ydiff)
    y=ya[n]+h/2*k1
    x=xa[n]+h/2
    k2=eval(ydiff)
    y=ya[n]+h/2*k2
    x=xa[n]+h/2
    k3=eval(ydiff)
    y=ya[n]+h*k3
    x=xa[n]+h
    k4=eval(ydiff)
    ya[n+1]=ya[n]+h/6*(k1+2*k2+2*k3+k4)
    s="     %d    |    %f  |    %f    |" %(n,xa[n],ya[n])
    print(s)

s="     %d    |    %f  |    %f    |" %(n+1,xa[n+1],ya[n+1])
print(s)
