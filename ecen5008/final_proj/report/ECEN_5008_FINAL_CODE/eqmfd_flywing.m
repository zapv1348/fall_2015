% build up eq mfd for flying wing
%
%   in trim, the velocity for (this) flying wing car (model) doesn't play a role
%
%
%   interesting, we find, empirically, that
%     beta and u1
%   look largely like
%     linear and quadratic
%   functions of
%     alat
%   up to at least 10 m/s^2  and reasonable after that
%   
% JH
% nov 15 boulder

disp_iter = 'off';
% disp_iter = 'iter';

v = 1;

% alat_s = (0:0.1:30)';
alat_s = (0:0.02:30)';
% alat_s = (0:0.1:20)';
% alat_s = (0:0.1:10)';
nn = length(alat_s);
nn2 = ceil(nn/2);

beta_s = 0*alat_s;
u1_s = beta_s;
Fy_s = beta_s;

u1_beta_0 = [ ];

for ii=1:length(alat_s)
  alat = alat_s(ii);
  omega = alat/v;
  [u1_beta, Fy] = trim_slcar([v omega], u1_beta_0, disp_iter);
  beta_s(ii) = u1_beta(2);
  u1_s(ii) = u1_beta(1);
  Fy_s(ii) = Fy;

  u1_beta_0 = u1_beta;
end

% save the equilibrium manifold
Alat_Beta_U1 = [ flipud([ -alat_s(2:end) -beta_s(2:end) u1_s(2:end) ]);
                 [ alat_s beta_s u1_s ] ];

save eqmfd_slc  Alat_Beta_U1

% beta vs alat looks fairly linear, and u1 vs alat quadratic
%
% do a least squares fit on alat=0 to 10
%
%   beta = beta_alat * alat
%   u1   = u1_alat2 * alat^2/2
alat = 10;
ii = min( find(alat_s >= alat) );
beta_alat = alat_s(1:ii)        \ beta_s(1:ii);
u1_alat2  = (alat_s(1:ii).^2/2) \ u1_s(1:ii);

save beta_u1_alat  beta_alat u1_alat2

figure
  plot(alat_s,Fy_s)
  grid on, zoom on
title('F_y vs a_l_a_t')

figure
  plot(-beta_s*180/pi,Fy_s)
  grid on,zoom on
title('F_y vs -\beta')

figure
  plot(alat_s,u1_s, alat_s, u1_alat2*alat_s.^2/2)
  grid on, zoom on
title('u_1 vs a_l_a_t')

figure
  plot(alat_s,beta_s*180/pi, alat_s, beta_alat*alat_s*180/pi)
  grid on, zoom on
title('\beta vs a_l_a_t')
