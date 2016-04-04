function rhs_xykk = c_rhs(xykk,X,Y,Iu,Iv,Iux,Iuy,Ivx,Ivy)
% Given x,y,k1,k2 in the vector xykk; calculate the RHS of the ODE system

x  = xykk(1); y  = xykk(2); k1 = xykk(3); k2 = xykk(4);
u  = Iu (y,x);                     % Interpolate to obtain u, v, etc. 
v  = Iv (y,x);                     % at the present location
ux = Iux(y,x);
uy = Iuy(y,x);
vx = Ivx(y,x);
vy = Ivy(y,x);  
g = 9.8 ;                          % With u,v,ux,uy,vx,vy known at the
k = sqrt(k1^2+k2^2);               % location x,y, we can calculate the
alpha = 0.5*sqrt(g/k^3);           % four components of the RHS   
rhs_xykk = [alpha*k1+u; alpha*k2+v; -(k1*ux+k2*vx); -(k1*uy+k2*vy)];  
end