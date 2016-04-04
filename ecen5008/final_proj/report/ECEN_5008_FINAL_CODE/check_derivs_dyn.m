% check derivatives of dynamics
% JH nov 15 boulder

% [dx, fxu_x,fxu_u, q_fxu_x_x,q_fxu_x_u,q_fxu_u_u] = dynamics_m(x,u,wt, q);
% [dx, A, B] = dynamics_m(x,u,wt);

% clear all, close all, clc

[NS, NI, NO, NW, NWL, NWC] = sys_sizes_m;

% flying wing

% x = [  v   gamma alpha omega]
xxx = [  10   0.1    0  0 ]';

% u = [ u1 u2 ]
uuu = [ 10  1  ]';

% wt = [ ]
wt   = [ ];
% wt = [ 1 2 ];   % test that the size check works


% x = [  v   gamma alpha omega]
zzz = [  0   0   0 0 ]';
zzz = rand(NS,1)/3

% u = [  u1 f_fx     f_rx    delta ]
% u = [ u1 u2 ]
vvv = [  0.5  0 ]';

eps_s = (-1:.001:1)';
% eps_s = (-.5:1e-4:-.4)';

fxu_s = [];
Dfxu_zv_s = [];

for ii = 1:length(eps_s)

  epsi = eps_s(ii);
%   disp(sprintf('eps : %g\n',epsi));


  % [dx, y, A, B] = dynamics_m(x,u,wt);
  [fxu, y, fxu_x, fxu_u] = dynamics_m(xxx+epsi*zzz, uuu+epsi*vvv, wt);

  fxu_s = [ fxu_s; fxu' ];

  Dfxu = [ fxu_x fxu_u ];

  Dfxu_zv_s = [ Dfxu_zv_s;  (Dfxu*[zzz; vvv])' ];

end

Dfxu_zv_s_fd = diff(fxu_s)  ./  repmat(diff(eps_s),1,length(fxu)) ;
Dfxu_zv_s_fd = [ Dfxu_zv_s_fd; Dfxu_zv_s_fd(end,:) ];

for ii =  1:length(xxx)

 if 0
  figure
    plot(eps_s, fxu_s(:,ii))
    grid on, zoom on
  title(sprintf('f_%d(x+\\epsilon z,u+\\epsilon v)',ii))
 end

  figure
    plot(eps_s, [ Dfxu_zv_s(:,ii) Dfxu_zv_s_fd(:,ii)])
    grid on, zoom on
  title(sprintf('d/d\\epsilon f_%d(x+\\epsilon z,u+\\epsilon v)',ii))

end



return
