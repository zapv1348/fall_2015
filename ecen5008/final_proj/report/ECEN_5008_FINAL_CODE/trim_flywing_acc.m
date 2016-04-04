function [alf_u1_u2, x, u, dx] = trim_flywing_acc(v_gam_vdot_gamdot,alf_u1_u2_0,display_iter)
% [alf_u1_u2, Fy] = trim_flywing(v_gamma, alf_u1_u2_0,display_iter)
%
%   trim the flying wing at a given v,gamma
%
%   display_iter: 'off', 'iter', 'final'
%     use 'off' silence iterations
%
% JH
% nov 15 boulder

if nargin < 3
  display_iter = 'iter';
end

opts = optimset('Display',display_iter,'TolX',1e-12,'TolFun',1e-10);

if nargin < 2
  alf_u1_u2_0 = [ ];
end

if isempty(alf_u1_u2_0)
  alf_u1_u2_0 = [ .1 1 0 ];
end

[alf_u1_u2, Fval, exitflag, out] = ...
    fsolve(@(alf_u1_u2)flywing_eq(alf_u1_u2, v_gamma), alf_u1_u2_0, opts);

if nargout > 1
  v =  v_gam_vdot_gamdot(1);
  gam = v_gam_vdot_gamdot(2);
  alf = alf_u1_u2(1);

  u1 = alf_u1_u2(2);
  u2 = alf_u1_u2(3);

  x = [v gam alf 0];
  u = [u1 u2];
  dx = dynamics_m([v gam alf 0],u);
end

if 1 ~= exitflag
  fprintf('exitflag: %d\n',exitflag);
  disp(out.message);
end

% keyboard


  % nested function
  function res = flywing_eq(alf_u1_u2,  v_gam_vdot_gamdot)

    v   = v_gam(1);
    gam = v_gam(2);
    alf = alf_u1_u2(1);
    th = gam + alf;

    u1 = alf_u1_u2(2);
    u2 = alf_u1_u2(3);

    dx = dynamics_m([v alf 0 th],[u1 u2]);

    res = dx(1:3);
  end

end
  