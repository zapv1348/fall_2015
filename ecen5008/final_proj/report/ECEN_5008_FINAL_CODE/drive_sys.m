% drive_sys - play with a system ... check proj operator, etc.
%
%   example:
%
%   play with the 
%
%     flying wing
%
%   NOTE: this script is a collection
%     of *many* little vignettes that
%     one might cut and paste into matlab
%     (many times with small changes)
%
%   John Hauser
%   apr12 ... mar 15
%   nov 15  boulder


% get system sizes
[NS, NI, NO, NW, NWL, NWC] = sys_sizes_m;

NS2 = NS*NS; NN12 = NS*(NS+1)/2;

% trim in forward flight
% disp_iter = 'iter';  % 'off';
% 
% v0 = 10;
% gam0 = 0;
% [alf_u1_u2,x0, u0] = trim_flywing([v0 gam0],[],disp_iter);
% 
% v1 = 10;
% gam1 = 10*pi/180;
% [~, x1, u1] = trim_flywing([v1 gam1],alf_u1_u2,disp_iter);
v0 = 10;
[~,x0,u0,~,~] = trim_acc([v0 0 0 0]); %%% trim around starting point
% linearize about trim
[~, ~, A0, B0] = dynamics_m(x0,u0);

% design an LTI feedback
Qreg = diag([10 50 20 40]);
Rreg = diag([.1 .1]);

[K,P,E] = lqr(A0,B0,Qreg,Rreg);
fprintf('E: %g %g %g %g\n',E);

Kr_ = [ K(1,:), K(2,:) ];

a  = 1;
dt = 0.05;
zmax = 60;
th2 = zmax*(2.5/16);
th1 = 0.1*th2;
th3 = 0.1*th2;
thalf = th1+th2+th3;
%section for half loop
T0 = 0:dt:th1;
T1 = (th1+dt):dt:(th1+th2);
T2 = (th1+th2+dt):dt:thalf;
gamma_t0 = zeros(1,length(T0));
gamma_t1 = (pi/th2)*(T1-th1);
gamma_t2 = pi*ones(1,length(T2));
Gamma = [gamma_t0 gamma_t1 gamma_t2];
T = [T0 T1 T2];
%T = T';
Gammadot=[zeros(1,length(T0)),pi/th2*ones(1,length(T1)),zeros(1,length(T2))];


tdes1 = 0.02;
tdes2 = 2.5;
tdes3 = 5.8;
tdes4 = 2.5;
tdes5 = 0.02;
tdes = tdes1+tdes2+tdes3+tdes4+tdes5;
desang = pi/4;
%section for descent
ta = tdes1+thalf;
tb = ta+tdes2;
tc = tb+tdes3;
td = tc+tdes4;
tf = td+tdes5-dt;

T0 = thalf+dt:dt:ta;
T1 = (ta+dt):dt:tb;
T2 = (tb+dt):dt:tc;
T3 = (tc+dt):dt:td;
T4 = (td+dt):dt:tf;

gamma1 = pi*ones(1,length(T0));
gamma2 = pi+((desang)/(tb-ta))*(T1-ta);
gamma3 = (pi+desang)*ones(1,length(T2));
gamma4 = (pi+desang)-(desang/(td-tc))*(T3-tc);
gamma5 = pi*ones(1,length(T4));

TS = [T0 T1 T2 T3 T4];
T=[T TS T+(tdes+thalf) TS+(tdes+thalf)];

Gamma1 = [gamma1 gamma2 gamma3 gamma4 gamma5];
Gamma=[Gamma Gamma1 flip(Gamma) (pi-Gamma1)];

Gammadot1 = [zeros(1,length(T0)) (desang/(tb-ta))*ones(1,length(T1)) zeros(1,length(T2)) -(desang/(td-tc))*ones(1,length(T3)) zeros(1,length(T4))]; 
Gammadot=[Gammadot, Gammadot1, -Gammadot  flip(Gammadot1)];

%%Gammadot = [gamma_t0 (pi/10)*ones(1,length(gamma_t1)) zeros(1,length(gamma_t2))];
V = v0*ones(1,length(T));
Vdot = zeros(1,length(T));

[Xdes,Udes,~] = build_Xdes_Udes(V,Gamma,Vdot,Gammadot,T);
Xdes = Xdes';
Udes=Udes';
Kr = repmat(Kr_,length(T),1);
T = T';
%%Xdes = repmat(x0,length(T0),1);
%%Udes = repmat(u0,length(T0),1);

Wt = [];
Wlt = [Xdes Udes];

% set up nonlin stuff

% we will use the 'sim' command (Simulink) to integrate ODEs (init value probs)
%   in order to use an old (obsolete) syntax, we need the following
load_system('LQR_KrS');
load_system('nonl_KS');
set_param('LQR_KrS',  'InitInArrayFormatMsg', 'None');
set_param('nonl_KS',  'InitInArrayFormatMsg', 'None');
save_system('LQR_KrS');
save_system('nonl_KS');


% these are some options for the integrator (sim)
opts = simset('Solver','ode45', ...
	      'AbsTol',1e-8,'RelTol',1e-6, ...
              'Trace','minstep');  % to catch blowup ... later (for Riccati)



% nonl_K implements the projection operator (+ cost evaluation)
%   note the (extended) input

del_x0 = [ -.1 -.02 0 0 ];
del_x0 = [ 0 0 0 0.1 ];
del_x0= [0 0 0 0];
x00 = x0' + del_x0;
ip_opts = simset(opts, 'InitialState', [x00 0]);
[Ti, aXi, aYi] = sim('nonl_KS', T, ip_opts, [T [Xdes Udes] Kr Wt Wlt]);
     %  the 'a' in aXi, aYi means 'augmented' (why?)

X1 = aXi(:,1:NS);
U1 = aYi(:,1:NI);
Y1 = aYi(:,NI+1:NI+NO);


figure,plot(T,[Xdes(:,1) X1(:,1)]),grid on,zoom on,title('V')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('Desired', 'Actual')

figure,plot(T,[Xdes(:,2) X1(:,2)]*180/pi),grid on, zoom on, title('\gamma')
xlabel('Time (sec)')
ylabel('Flight Path Angle (Deg)')
legend('Desired', 'Actual')

figure,plot(T,[Xdes(:,3) X1(:,3)]*180/pi),grid on, zoom on, title('\alpha')
xlabel('Time (sec)')
ylabel('Angle of Attack (Deg)')
legend('Desired', 'Actual')

figure,plot(T,[Xdes(:,4) X1(:,4)]*180/pi),grid on, zoom on, title('\omega')
xlabel('Time (sec)')
ylabel('Omega (Deg/sec)')
legend('Desired', 'Actual')



figure,plot(T,[Udes(:,1) U1(:,1)]),grid on,zoom on,title('u_1')
xlabel('Time (sec)')
ylabel('u1 (N)')
legend('Desired', 'Actual')

figure,plot(T,[Udes(:,2) U1(:,2)]),grid on,zoom on,title('u_2')
xlabel('Time (sec)')
ylabel('u2 N')
legend('Desired', 'Actual')



% save some stuff for traj optimization

Tdes = T;

X0 = X1;
U0 = U1;

% make sure Qreg and Rreg get updated in  QR_params.h  !!

save des_traj  Tdes Xdes Udes  T X0 U0
v = @(t) interp1(T,X1(:,1),t);
gamma = @(t) interp1(T,X1(:,2),t);


%%% Integrate to get x and z
vdes = @(t)interp1(T,Xdes(:,1),t);
gammades = @(t) interp1(T,Xdes(:,2),t);
kinemat = @(t,x) [v(t)*cos(gamma(t)); -v(t)*sin(gamma(t))];
kinematdes = @(t,x) [vdes(t)*cos(gammades(t)); -vdes(t)*sin(gammades(t))];
[Ti,XZ] = ode45(kinemat,T,[0 0],odeset('AbsTol',1e-8,'RelTol',1e-6));
[Tides, XZdes] = ode45(kinematdes,T,[0 0],odeset('AbsTol',1e-8,'RelTol',1e-6));
figure,plot(XZdes(:,1),-XZdes(:,2),XZ(:,1),-XZ(:,2)),grid on,zoom on,axis equal,title('Wing Path')
xlabel('X  (m)');
ylabel('-Z (m)');
legend('Desired','Actual')


return

