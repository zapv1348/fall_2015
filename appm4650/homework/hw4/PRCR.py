import numpy as np

#predictor:
#y_{i+1}=y_i+\frac{h}{24}(55f_i-59*f_{i-1}+37f_{i-2}-9f_{i-3})

#corrector:
#y_{i+1}=y_i+\frac{h}{24}(9f_{i+1}+19f_i-5f_{i-1}+f_{i-2})

ydiff="x+y"

seed="x**2/2"

x0=0.0
y0=0.0

xmax=0.5
h=0.1

def my_range(start,end,step):
    while start<=end:
        yield start
        start+=step




#start with step of 0.025
x=x0
y=y0
f0=eval(ydiff)

x=0.025
y=eval(seed)
f1=eval(ydiff)

x=0.05
y=eval(seed)
f2=eval(ydiff)

x=0.075
y=eval(seed)
f3=eval(ydiff)

y1=y+0.025/24*(55*f3-59*f2+37*f1-9*f0)

#now we have the first predictor value, let us correct with h=0.033
x=1.0/30.0
y=eval(seed)
f1=eval(ydiff)

x=2.0/30.0
y=eval(seed)
ye=y
f2=eval(ydiff)

x=0.1
y=y1
f3=eval(ydiff)

y1=ye+0.033/24*(9*f3+19*f2-5*f1+f0)

#now we do the second predictor with h=0.05
#thus, we need f values at x=0,.05,0.1,0.15


x=0.05
y=eval(seed)
f2=eval(ydiff)

x=0.1
y=y1
f4=eval(ydiff)

x=0.15
y=eval(seed)
f5=eval(ydiff)

y2=y+h/24*(55*f5-59*f4+37*f2-9*f0)

#still using h=0.05
x=0.2
y=y2
f6=eval(ydiff)

y2=0.15+0.05/25*(9*f6+19*f5-5*f4+f2)

y=y2
f6=eval(ydiff)

x=0.1
y=y1
f8=eval(ydiff)

x=0.25
y=eval(seed)
f7=eval(ydiff)

y3=y+.05/24*(55*f7-59*f6+37*f5-9*f4)

y3=y2+0.1/24*(9*y3+19*f6-5*f8+0)
x=0.3
y=y3
f9=eval(ydiff)

y4=y3+0.1/24*(55*f9-59*f6+37*f8-0)
x=0.4
y=y4
f10=eval(ydiff)

y4=y3+0.1/24*(9*f10+19*f9-5*f6+f8)
y=y4
f10=eval(ydiff)

y5=y4+0.1/24*(55*f10-59*f9+37*f6-9*f8)
x=0.5
y=y5
f11=eval(ydiff)

y5=y4+0.1/24*(9*f11+19*f10-5*f9+f6)

x=np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
y=np.array([0, y1, y2, y3, y4, y5])

print("iteration|       x     |       y        |")

for n in my_range(0,5,1):
    s="    %d    |  %f   |    %f    |" %(n,x[n],y[n])
    print(s)
