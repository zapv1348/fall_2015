import numpy as np

ydiff="y+x"

y0=0.0
x0=0.0

xmax=0.5
h=0.1

yest=np.array([y0,0,0,0,0,0])

#improved euler method is:
#y_{j+1}=y_j+\frac{h}{2}(f(x_j,y_j)+f(x_{j+1},y_j+hf(x_j,y_j)))


def my_range(start,end,step):
    while start<=end:
        yield start
        start+=step

y=yest[0]


r="iteration|       x        |       y        |"
print(r)

for x in my_range(x0,xmax-h,h):
    it=x*10
    m1=eval(ydiff)
    x=x+h
    y=y+m1*h
    m2=eval(ydiff)
    yest[it+1]=yest[it]+h*(m1+m2)/2
    s="    %d    |    %f    |    %f    |" %(it,x-h,yest[it])
    print(s)

s="    %d    |    %f    |    %f    |" %(it+1,x,yest[it+1])
print(s)

