function [ v2prime ] = G2( u1, v1, v2, Pr )
%Differential equation for v2( v2 = G' from original problem)
v2prime = -Pr * u1 * v2 /2;

end

