% check the second derivatives of dynamics
% FB nov 11 boulder
% JC jul 15 boulder
% JH nov 15 boulder

% [dx, fxu_x,fxu_u, q_fxu_x_x,q_fxu_x_u,q_fxu_u_u] = dynamics_m(x,u,wt, q);
% [dx, A, B] = dynamics_m(x,u,wt);

% clear all, close all, clc

[NS, NI, NO, NW, NWL, NWC] = sys_sizes_m;

% look at some equilibrium states
v = 10; alat = 5;
omega = alat/v;
[u1_beta,Fy] = trim_slcar([v omega],[],'off');
beta = u1_beta(2);
u1 = u1_beta(1);
u2 = 0;

% sliding car

% x = [  v   beta omega ]
% xxx = [  10    0    0   ]';
xxx = [  v   beta omega ]';

% u = [ u1 u2 ]
% uuu = [ 0  0  ]';
uuu = [ u1 u2 ]';

% wt = [ ]
wt   = [ ];
% wt = [ 1 2 ];   % test that the size check works

% x = [  v   beta omega ]
zzz = [  0   0.1   0    ]';
zzz = rand(NS,1)/NS;

% u = [  u1 f_fx     f_rx    delta ]
% u = [ u1 u2 ]
vvv = [  1  0 ]';
vvv = rand(NI,1)/NI;


eps_s = (-1:.001:1)';
% eps_s = (-.5:1e-4:-.4)';

%   [ v  beta omega ]
q = [ 0  1   0 ]';
q = rand(NS,1)/NS;

gxu_s = [];
Dgxu_zv_s = [];

for ii = 1:length(eps_s)

  epsi = eps_s(ii);
%   disp(sprintf('eps : %g\n',epsi));

  % [~, ~, A, B, Q,S,R] = dynamics_m(x,u,wt,q);
  [~, ~, fxu_x, fxu_u, gx_x, gx_u, gu_u] = ...
                               dynamics_m(xxx+epsi*zzz, uuu+epsi*vvv, wt, q);

  gxu = (q'*[fxu_x fxu_u] )';

  gxu_s = [ gxu_s; gxu' ];

  Dgxu = [gx_x gx_u; gx_u' gu_u];

  Dgxu_zv_s = [ Dgxu_zv_s;  (Dgxu*[zzz; vvv])' ];

end

Dgxu_zv_s_fd = diff(gxu_s)  ./  repmat(diff(eps_s),1,length(gxu)) ;
Dgxu_zv_s_fd = [ Dgxu_zv_s_fd; Dgxu_zv_s_fd(end,:) ];

for ii =  1:length(gxu)

 if 0
  figure
    plot(eps_s, gxu_s(:,ii))
    grid on, zoom on
  title(sprintf('g_%d(x+\\epsilon z,u+\\epsilon v)',ii))
 end

  figure
    plot(eps_s, [ Dgxu_zv_s(:,ii) Dgxu_zv_s_fd(:,ii)])
    grid on, zoom on
  title(sprintf('d/d\\epsilon g_%d(x+\\epsilon z,u+\\epsilon v)',ii))

end



return

