function [ u3prime ] = F3( u1, u2, u3 )
%Differential equation for u3, (u3 = F'' from original problem)
u3prime = -.5 * u1 * u3;


end

