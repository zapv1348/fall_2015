#numpy is what we will use for the majority of our calculations
import numpy as np

#library for ploting
import matplotlib.pyplot as plt

#library for exiting if stuff is bunk
import sys


#delta for the fizzle
delta=0.2

#diff eq to approximate with rk4, x is theta
dfizzle="0.2*np.exp(x)-x"


#initial condition
x0=0.0
sig0=0.0

#step size to use
h=0.1

#steps
step=200

#a nice thing for generating for loops
def my_range(start,end,step):
    while start<=end:
        yield start
        start+=step

#initialize arrays
omega=np.zeros((step,1))
sig=np.zeros((step,1))

omega[0]=x0
sig[0]=sig0


#now to do rk4 for this problem
for i in my_range(0,step-2,1):
    x=omega[i]
    sig[i+1]=sig[i]+h
    k1=eval(dfizzle)
    x=omega[i]+h*k1/2
    k2=eval(dfizzle)
    x=omega[i]+h*k2/2
    k3=eval(dfizzle)
    x=omega[i]+h*k3
    k4=eval(dfizzle)
    omega[i+1]=omega[i]+h/6*(k1+2*k2+2*k3+k4)

#print(omega[step-2])
#print(omega[step-1])


#gotta install the plotting stuff first
#plt.plot(sig,omega)
#plt.show()


####################################
#Now we will do the root finding method to find the theta that
#gives the approximate llong term solution

#Brent's method cause doing new things is good for you

approx="np.exp(x)/x-5"

a=0.2
b=2.0



x=a
fa=eval(approx)
x=b
fb=eval(approx)


if fa*fb>=0:
    sys.exit()

fa1=np.absolute(fa)
fb1=np.absolute(fb)
if fa1<fb1:
    r=a
    a=b
    b=r
    r=fa
    fa=fb
    fb=r
c=a
fc=fa

s=0
fs=50

delt=0.000001

flag=1
d=0

while (np.absolute(fb)>=0.0000001) & (np.absolute(fs)>=0.0000001):
    if ((fa!=fc) & (fb!=fc)):
        s=a*fb*fc/((fa-fb)*(fa-fc))+b*fa*fc/((fb-fa)*(fb-fc))+c*fa*fb/((fc-fa)*(fc-fb))
    else:
        s=b-fb*(b-a)/(fb-fa)
    if  (~((3*a+b)/4<s<b))|(flag&(np.absolute(s-b)>=(np.absolute(b-c)/2)))|(~flag&(np.absolute(s-b)>=np.absolute(c-d)/2))|(flag&(np.absolute(b-c)<delt))|(flag&(np.absolute(c-d)<delt)):
        s=(a+b)/2
        flag=1
    else:
        flag=0
    x=s
    fs=eval(approx)
    d=c
    c=b
    x=c
    fc=eval(approx)
    if (fa*fs)<0:
        b=s
        x=b
        fb=eval(approx)
    else:
        a=s
        x=a
        fa=eval(approx)
    if np.absolute(fa)<np.absolute(fb):
        r=a
        a=b
        b=r
        r=fa
        fa=fb
        fb=r

if fb<fs:
    print (b)
else:
    print (s)


#########################
#here, we will be plotting the early and late solutions with the numeric solution.
shortt="(-.25)*(np.exp(-.8*y)-1)"
step=30
early=np.zeros((step,1))
yx=np.zeros((step,1))
for n in my_range(0,step-2,1):
    yx[n+1]=yx[n]+h
    y=yx[n]
    early[n]=eval(shortt)

y=yx[n+1]
early[n+1]=eval(shortt)

step=135
late=np.ones((step,1))*b
yr=np.linspace(6.5,20,step)



line1=plt.plot(sig,omega,label="Whole",linewidth=1.5,linestyle="--")
line2=plt.plot(yx,early,label="Early",linewidth=1.5)
line3=plt.plot(yr,late,label="Late",linewidth=1.5)
plt.title("Early and Late Solutions sigma vs omega",fontsize=24)
plt.xlabel("sigma",fontsize=24)
plt.ylabel("theta",fontsize=24)
plt.legend()
plt.show()
