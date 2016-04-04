function [ eta ] = Interpolate( x0, x1, x2, f0, f1, f2 , xi)
%Takes in three points (x0, f0), (x1, f1), (x2, f2) and uses Lagrange
% polynomials to approximate the function value at the location xi
eta = f0 * (xi - x1) * (xi - x2)/((x0 - x1) * (x0 - x2)) + f1 * (xi - x0) * (xi - x2)/((x1 - x0) * (x1 - x2)) + f2 * (xi - x0) * (xi - x1)/((x2 - x0) * (x2 - x1));

end

