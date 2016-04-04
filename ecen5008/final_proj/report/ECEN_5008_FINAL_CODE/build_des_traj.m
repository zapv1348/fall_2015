% specify desired
%
%     v(t)
%
%   and
%
%     alat(t)
%
%   (+ derivs)
%
%   and then, using the equilibrium manifold,
%     build up a desired (x,u)(.) [almost] trajectory
%
% JH
% nov 15 boulder

% to defined a 'desired traj',
%   we will take a 
%
%     quasi-static approach
%
%   wherein we imagine (falsely) that, at each point in time,
%   we are on a constant speed circle determined by the
%
%     current v and alat
%

% here's we will use the (pre-computed) eq mfd to build a desired trajectory
%   recall that beta and u1 do *not* depend on v

load eqmfd_slc  % get Alat_Beta_U1


% our 'play example' will use a
%   constant v  and a
%   tanh     alat
% curve

Ttrans = 4;  % time to transition from alat=0 to alat=alat_max using  tanh
Ttrans = 2;  % time to transition from alat=0 to alat=alat_max using  tanh
Ttrans = 1;  % time to transition from alat=0 to alat=alat_max using  tanh

t1 = 3*Ttrans;

% set up a function space: choose a discretation level & (init &) final time
dt = 0.01;
T0 = (-t1:dt:t1)';

% choose v, alat
v0 = 10;
alat_max = 10;
alat_max = 20;



% tanh 4 = 0.99933
% tanh 3 = 0.99505
Alat_des  = alat_max*tanh(4/Ttrans * T0);
dAlat_des = alat_max*(1 - tanh(4/Ttrans * T0).^2)*4/Ttrans;

Beta_des = interp1(Alat_Beta_U1(:,1),Alat_Beta_U1(:,2), Alat_des);
U1_des   = interp1(Alat_Beta_U1(:,1),Alat_Beta_U1(:,3), Alat_des);

% use forward difference approximation
dBeta_des = diff(Beta_des) ./ diff(T0);
dBeta_des = [ dBeta_des; dBeta_des(end) ];

% v is constant (for now)
V_des = repmat(v0,size(T0));
dV_des = 0*V_des;

Omega_des  = Alat_des ./ V_des - dBeta_des;
% Omega_des  = Alat_des ./ V_des;  % be quasi-static ...
% use forward difference approximation
dOmega_des = diff(Omega_des) ./ diff(T0);
dOmega_des = [ dOmega_des; dOmega_des(end) ];


U2_des = dOmega_des;
% U2_des = 0*dOmega_des;

Xdes = [ V_des Beta_des Omega_des ];
Udes = [ U1_des U2_des ];

Tdes = T0;

% in principle, we'd write a dedicated script for this
save des_traj  Tdes Xdes Udes Alat_des
