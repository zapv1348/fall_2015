close all
clear all
t=linspace(0,2*pi,1000);

x1=cos(t);
x2=sin(t);

x=[x1;x2];

A=[1 3;0 1];
B=[2 -4; 1 4];
plot(x1,x2)
hold on
A1=A*x;
plot(A1(1,:),A1(2,:));
B1=B*x;
hold on
plot(B1(1,:),B1(2,:));

z=linspace(1,1,1000);
r=linspace(1,-1,1000);
x1=[z r -z -r];
x2=[-r z r -z];
x=[x1;x2];
figure
plot(x1,x2)
A2=A*x;
hold on
plot(A2(1,:),A2(2,:));
B2=B*x;
hold on
plot(B2(1,:),B2(2,:));