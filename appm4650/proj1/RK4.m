function [ sigmanew ] = RK4( theta, sigma, h )
%Run RK4 runs a single iteration
%   Detailed explanation goes here
k1 = h * F(theta);
k2 = h * F(theta + h/2);
k3 = h * F(theta + h/2);
k4 = h * F(theta + h);
sigmanew = sigma + k1/6 + k2/3 + k3/3 + k4/6;

end

