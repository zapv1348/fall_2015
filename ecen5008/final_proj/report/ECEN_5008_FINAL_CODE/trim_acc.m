function [alf_u1_u2,x,u,dx,error] =  trim_acc(v_gam_vdot_gamdot)
%constants for ereyting
p = 1.2;
S = 0.61;
Cd0 =0.1716;
Cd2 = 2.395;
Cl1 = 3.256;
c_bar = 0.5;
Cm1 = -0.1;
J = 0.24;
lz = 0.31;
o_m = 1.0/12.0;
g = 0.6;
%initial guess for fsolve
alf_u1_u2_0 = [ .1 1 0 ];

%array of paramaters
params = [p S Cd0 Cd2 Cl1 c_bar Cm1 J lz o_m g];
%This makes a function
fun1 = @(alf_u1_u2)acc_trim_dyn(v_gam_vdot_gamdot,alf_u1_u2,params);
opts = optimset('TolX',1e-12,'TolFun',1e-10);
%solve our dynamics for  where they are = to 0
alf_u1_u2 = fsolve(fun1,alf_u1_u2_0, opts);
%organize it into vectors to return
x = [v_gam_vdot_gamdot(1) v_gam_vdot_gamdot(2) alf_u1_u2(1) v_gam_vdot_gamdot(4)]';
u = [alf_u1_u2(2) alf_u1_u2(3)]';
dx = dynamics_m(x,u);
error = [v_gam_vdot_gamdot(3) v_gam_vdot_gamdot(4) 0 0]'-dx;
end