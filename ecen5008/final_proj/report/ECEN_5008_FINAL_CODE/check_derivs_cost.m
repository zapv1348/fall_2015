% check derivatives of cost function
%
%   flight path cfv system
%
% uses flight path frame, parametrized by (chi, gamma)


% [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(x, u, wlt)

% clear all, close all, clc

pronto_paths

% x = [  v  chi  gam  ]
xxx = [  10  0   .2  ]';

% u = [  ax ay az ]
uuu = [  0  0 -9.81 ]';

% wlt   = [0 0 0 0 0 0 ];
wlt = [ xxx; uuu ];

% x = [  v  chi  gam  ]
zzz = 0*[  1  .1  .1 ]';
% zzz = rand(4,1)/3

% u = [  f_fx     f_rx    delta ]
vvv = [ 0 1 0]';

eps_s = (-1:.001:1)';
% eps_s = (-.5:1e-4:-.4)';

lxu_s = [];
Dlxu_s = [];
Dlxu_zv_s = [];
D2lxu_zv_s = [];
for ii = 1:length(eps_s)

  epsi = eps_s(ii);
%   disp(sprintf('eps : %g\n',epsi));


% [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(x, u, wlt)
  [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(xxx+epsi*zzz,...
                                                          uuu+epsi*vvv, wlt);

  lxu_s = [ lxu_s; lxu ];

  Dlxu = [ lxu_x lxu_u ];
  Dlxu_s = [Dlxu_s; Dlxu];
  Dlxu_zv_s = [ Dlxu_zv_s;  Dlxu*[zzz; vvv] ];

  % Build hessian 
  D2lxu = [ lxu_x_x lxu_x_u; 
            lxu_x_u' lxu_u_u];

  D2lxu_zv_s = [D2lxu_zv_s; (D2lxu*[zzz;vvv])' ];


end

Dlxu_zv_s_fd = diff(lxu_s)  ./  repmat(diff(eps_s),1,length(lxu)) ;
Dlxu_zv_s_fd = [ Dlxu_zv_s_fd; Dlxu_zv_s_fd(end,:) ];

D2lxu_zv_s_fd = diff(Dlxu_s) ./ repmat(diff(eps_s),1,length(Dlxu)) ;
D2lxu_zv_s_fd = [D2lxu_zv_s_fd; D2lxu_zv_s_fd(end,:) ];
%plot first derivatives

figure
  plot(eps_s, [ Dlxu_zv_s Dlxu_zv_s_fd ])
  grid on, zoom on
title(sprintf('d/d\\epsilon l(x+\\epsilon z,u+\\epsilon v)'))


for ii =  1:length(xxx)

  figure
    plot(eps_s, [ D2lxu_zv_s(:,ii) D2lxu_zv_s_fd(:,ii)])
    grid on, zoom on
  title(sprintf('d/d\\epsilon dl_{x%d}(x+\\epsilon z,u+\\epsilon v)',ii))

end

for ii =  1:length(uuu)

  figure
    plot(eps_s, [ D2lxu_zv_s(:,ii+3) D2lxu_zv_s_fd(:,ii+3)])
    grid on, zoom on
  title(sprintf('d/d\\epsilon dl_{u%d}(x+\\epsilon z,u+\\epsilon v)',ii))

end



return
