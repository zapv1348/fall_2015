function [ I ] = Simpson( xi,h, N)
%Apply Simpson's method to F(x) [F is the differential equation defined in
%F.m]
%   Simpson's method is used to integrate the differential equation to find
%   sigma explosion.
sum = 0;
for i = 1:N/2
    sum = sum + F(xi + h * (2*i -2)) + 4 *F(xi + h*(2*i -1)) + F(xi + h *(2*i));
    
end
I = sum * h/3;

