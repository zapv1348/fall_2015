function [ y ] = LT( th )
%After finding the asymptotic value from integration, it is used for the 
%long term solution, which takes in a vector and returns a vector with
% each element pased through the long term approximation
for i = 1:length(th)
    y(i) = 1.3591 - 1/(exp(th(i)));
end

