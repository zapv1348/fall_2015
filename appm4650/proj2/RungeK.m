function [ m1 ] = RungeK( u1i, u2i, u3i,v1i, v2i, h, N, Pr )
%Run RK4 using Rk4 helper function for N iterations returns a matrix
% where the columns are eta, F, F', F'', G, G'. For each value of Pr, this
% function was called (with guesses for F''(0) and G'(0)) this matrix was 
% then used to generate the plots and to interpolate the values of etam and
% etat
%
m1 = zeros (N, 6);
m1(1,1) = 0;
m1(1,2) = u1i;
m1(1,3) = u2i;
m1(1,4) = u3i;
m1(1,5) = v1i;
m1(1,6) = v2i;
eta = 0;
u1 = u1i;
u2 = u2i;
u3 = u3i;
v1 = v1i;
v2 = v2i;
for i = 2:N
    uv = RK4(u1, u2, u3, v1, v2, h, Pr);
    u1 = uv(1);
    u2 = uv(2);
    u3 = uv(3);
    v1 = uv(4);
    v2 = uv(5);
    eta = eta + h;
    m1(i, 1) = eta;
    m1(i, 2) = u1;
    m1(i, 3) = u2;
    m1(i, 4) = u3;
    m1(i, 5) = v1;
    m1(i, 6) = v2;




end

