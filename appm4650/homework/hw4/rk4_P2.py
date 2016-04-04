import numpy as np

u1diff="3*u1+2*u2-(2*t**2+1)*np.exp(2*t)"

u2diff="4*u1+u2+(t**2+2*t-4)*np.exp(2*t)"

u1act="1/3*np.exp(5*t)-1/3*np.exp(-t)+np.exp(2*t)"
u2act="1/3*np.exp(5*t)+2/3*np.exp(-t)+t**2*np.exp(2*t)"

u10=1.0
u20=1.0
t0=0.0

tmax=1.0
h=0.2

u1a=np.array([u10, 0, 0, 0, 0, 0])
u2a=np.array([u20, 0, 0, 0, 0, 0])
ta=np.array([t0, t0+h, t0+2*h, 3*h, 4*h, 5*h])


u1e=np.array([0.0,0,0,0,0,0])
u2e=np.array([0.0,0,0,0,0,0])


def my_range(start,end,step):
    while start<=end:
        yield start
        start+=step

for t in my_range(t0,tmax-h,h):
    n=t*5
    u1e[n]=eval(u1act)
    u2e[n]=eval(u2act)
    u1=u1a[n]
    u2=u2a[n]
    t=ta[n]
    k1u1=eval(u1diff)
    k1u2=eval(u2diff)
    u1=u1a[n]+h/2*k1u1
    u2=u2a[n]+h/2*k1u2
    t=ta[n]+h/2
    k2u1=eval(u1diff)
    k2u2=eval(u2diff)
    u1=u1a[n]+h/2*k2u1
    u2=u2a[n]+h/2*k2u2
    t=ta[n]+h/2
    k3u1=eval(u1diff)
    k3u2=eval(u2diff)
    u1=u1a[n]+h*k3u1
    u2=u1a[n]+h*k3u2
    t=ta[n]+h
    k4u1=eval(u1diff)
    k4u2=eval(u2diff)
    u1a[n+1]=u1a[n]+h/6*(k1u1+2*k2u1+2*k3u1+k4u1)
    u2a[n+1]=u1a[n]+h/6*(k1u2+2*k2u2+2*k3u2+k4u2)

t=1
u1e[5]=eval(u1act)
u2e[5]=eval(u2act)

r="iteration |       t      | u1 estimate  |  u1 actual  | u2 estimate  |   u2 actual"
print(r)
for i in my_range(0,5,1):
    s="     %d    |    %f  |   %f   |  %f   |   %f   |   %f   " %(i+1,ta[i],u1a[i], u1e[i],u2a[i], u2e[i])
    print(s)
