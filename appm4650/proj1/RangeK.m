function [ m1 ] = RangeK( thetai, sigmai, h, N )
%Uses other RK4 function as a helper to run RK4 the specified number of iterations
% This code was used to create a 2 column matrix which was plotted as the
% numeric solution
m1 = zeros(N, 2);
m1(1,1) = thetai;
m1(1,2) = sigmai;
th = thetai;
sig = sigmai;
for i = 2:N
    sig = RK4(th, sig, h);
    th = th + h;
    m1(i, 1)= th;
    m1(i, 2) = sig;
    

    

end

