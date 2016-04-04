function F = acc_trim_dyn(v_gam_vdot_gamdot,alf_u1_u2,params)
    %
    v2 = v_gam_vdot_gamdot(1)^2;
    %parameter functions
    D = 0.5*params(1)*params(2)*v2*(params(3)+params(4)*alf_u1_u2(1)^2);
    L = 0.5*params(1)*params(2)*v2*params(5)*alf_u1_u2(1);
    M = 0.5*params(1)*params(2)*params(6)*v2*params(7)*alf_u1_u2(1);
    %the dynamic equations setup so that we can solve for F=0
    F(1) = v_gam_vdot_gamdot(3) + D*params(10)+params(11)*sin(v_gam_vdot_gamdot(2))-params(10)*(alf_u1_u2(2)*cos(alf_u1_u2(1))-alf_u1_u2(3)*sin(alf_u1_u2(1)));
    F(2) = -v_gam_vdot_gamdot(1)*v_gam_vdot_gamdot(4)+params(10)*L - params(11)*cos(v_gam_vdot_gamdot(2))+params(10)*(alf_u1_u2(2)*sin(alf_u1_u2(1))-alf_u1_u2(3)*cos(alf_u1_u2(1)));
    F(3) = M/(params(8))+(params(9)/params(8))*alf_u1_u2(3);
end