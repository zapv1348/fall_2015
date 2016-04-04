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

v0 = 10;
gam0 = 0;
[alf_u1_u2, x0, u0] = trim_flywing([v0 gam0]);

% linearize about trim
[~, ~, A0, B0] = dynamics_m(x0,u0);

% design an LTI feedback
Qreg = diag([20 40 1 40]);
Rreg = diag([1 1]);

[K,P,E] = lqr(A0,B0,Qreg,Rreg);
E

Kr_ = [ K(1,:), K(2,:) ];

t1 = 10;
T0 = (0:.01:10)';

Kr = repmat(Kr_,length(T0),1);

Xdes = repmat(x0,length(T0),1);
Udes = repmat(u0,length(T0),1);

Wt = [];
Wlt = [ Xdes Udes ];

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
x00 = x0 + del_x0;
ip_opts = simset(opts, 'InitialState', [x00 0]);
[Ti, aXi, aYi] = sim('nonl_KS', T0, ip_opts, [T0 [Xdes Udes] Kr Wt Wlt]);
     %  the 'a' in aXi, aYi means 'augmented' (why?)

X1 = aXi(:,1:NS);
U1 = aYi(:,1:NI);
Y1 = aYi(:,NI+1:NI+NO);


figure,plot(T0,[Xdes(:,1) X1(:,1)]),grid on,zoom on,title('V')
figure,plot(T0,[Xdes(:,4)-Xdes(:,2) X1(:,4)-X1(:,2)]*180/pi),grid on,zoom on,title('\gamma')

if 0
figure,plot(T0,[Xdes(:,2) X1(:,2)]*180/pi),grid on,zoom on,title('\alpha')
figure,plot(T0,[Xdes(:,3) X1(:,3)]*180/pi),grid on,zoom on,title('\omega')
figure,plot(T0,[Xdes(:,4) X1(:,4)]*180/pi),grid on,zoom on,title('\theta')
end

figure,plot(T0,[Udes(:,1) U1(:,1)]),grid on,zoom on,title('u_1')
figure,plot(T0,[Udes(:,2) U1(:,2)]),grid on,zoom on,title('u_2')



return

% from eqmfd_slcar.m
%   we saw that
%     beta(alat) is approx linear:     beta_alat * alat
%     u1(alat)   is approx quadratic:  u1_alat2  * alat^2/2
%   load those coefficients
load  beta_u1_alat  % get beta_alat u1_alat2


% specify desired
%
%     v(t) (& derivative)
%
%   and
%
%     alat(t)
%
%   and then use the simplified beta/u1 model(s) to
%     build up a desired (x,u)(.) [almost] trajectory
%

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

% our 'play example' will use a
%   constant v  and a
%   tanh     alat
% curve

Ttrans = 2;  % time to transition from alat=0 to alat=alat_max using  tanh
Ttrans = 4;  % time to transition from alat=0 to alat=alat_max using  tanh

t1 = 3*Ttrans;

% set up a function space: choose a discretation level & (init &) final time
dt = 0.01;
T0 = (-t1:dt:t1)';


alat_max = 10;

% tanh 4 = 0.99933
% tanh 3 = 0.99505
Alat_des  = alat_max*tanh(4/Ttrans * T0);
dAlat_des = alat_max*(1 - tanh(4/Ttrans * T0).^2)*4/Ttrans;


v0 = 10;
V_des = repmat(v0,size(T0));
dV_des = 0*V_des;

% our empirical linear model
Beta_des  = beta_alat * Alat_des;
dBeta_des = beta_alat * dAlat_des;


Omega_des  = Alat_des ./ V_des - dBeta_des;
% use forward difference approximation
dOmega_des = diff(Omega_des) ./ diff(T0);
dOmega_des = [ dOmega_des; dOmega_des(end) ];

U1_des = u1_alat2 * Alat_des.^2/2;   %maybe: + dV_des./cos(Beta_des);

U2_des = dOmega_des;

Xdes = [ V_des Beta_des Omega_des ];
Udes = [ U1_des U2_des ];

Tdes = T0;

% in principle, we'd write a dedicated script for this
save des_traj  Tdes Xdes Udes


x0 = [ v0 0 0 ];
u0 = [ 0 0 ];


% let's get the linearization at (x0,u0)
%
%   use:  doc dynamics_m
%     to see the calling sequence
[ ~, ~, A0, B0 ] = dynamics_m(x0,u0);


% grab the parameters we set in
%   QR_params.h
% for use with a time-varying LQR controller
%
%   use  doc QR_params_m  to see syntax
[Qr, Rr] = QR_params_m;

if 1
  % or PLAY with LQR design !!  ... be 
  Qr = diag([ 10 100 10 ]);
  Rr = diag([ 1 1 ]);
end

% for now, try a time invariant K
[K,~,E] = lqr(A0, B0, Qr, Rr);

% since NI > 1, we need to flatten K
K_ = [];
for ii=1:size(K,1)
  K_ = [ K_, K(ii,:) ];
end

Kr = repmat(K_,length(T0),1);    % Kr(t) for nonl sim


% set the 'desired traj' for tracking (as opposed to for the cost)
Xi = Xdes;
Ui = Udes;

% used in dynamics and cost
Wt = [ ];
Wlt = [ Xdes Udes ];


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

x00 = Xdes(1,:);
ip_opts = simset(opts, 'InitialState', [x00 0]);
[Ti, aXi, aYi] = sim('nonl_KS', T0, ip_opts, [T0 [Xi Ui] Kr Wt Wlt]);
     %  the 'a' in aXi, aYi means 'augmented' (why?)

X1 = aXi(:,1:NS);
U1 = aYi(:,1:NI);
Y1 = aYi(:,NI+1:NI+NO);


figure,plot(T0,[Xdes(:,1) X1(:,1)]),grid on,zoom on,title('V')
figure,plot(T0,[Xdes(:,2) X1(:,2)]*180/pi),grid on,zoom on,title('\beta')
figure,plot(T0,[Xdes(:,3) X1(:,3)]*180/pi),grid on,zoom on,title('\omega')

figure,plot(T0,[Alat_des Y1(:,2)]),grid on,zoom on,title('a_l_a_t')

figure,plot(T0,[U1_des U1(:,1)]),grid on,zoom on,title('u_1')
figure,plot(T0,[U2_des U1(:,2)]),grid on,zoom on,title('u_2')

% use anonymous functions with ode45 to integrate the kinematics
v = @(t) interp1(T0,X1(:,1),t);
beta = @(t) interp1(T0,X1(:,2),t);
omega = @(t) interp1(T0,X1(:,3),t);
kinemat = @(t,x) [ v(t)*cos(x(3)+beta(t)); v(t)*sin(x(3)+beta(t)); omega(t) ];

psi0 = -65*pi/180;
[Ti,XYPsi] = ode45(kinemat,T0,[0 0 psi0],odeset('AbsTol',1e-8,'RelTol',1e-6));

figure,plot(XYPsi(:,2),XYPsi(:,1)),grid on,zoom on,axis equal,title('ground track')

return




% here's some old stuff ... might be some pearls ...


% build up some kind of a "desired trajectory"
%
%   here, I've chosen an analytic expression,
%   but there are many other kinds of things we might do.
%
%   eventually, the desired traj can be used in L2 optimization
%   here, I'm only trying to get a curve to *play* around with
%   (like in class)

% define a desired phi traj

tanhT2 = tanh(T0-t1/2);
Phi_0   = tanhT2.*(1-tanhT2.^2);
Phi_0_max = max(Phi_0);

amp = ( 20 *pi/180 );

Phi_des   = amp*(tanhT2.*(1-tanhT2.^2))/Phi_0_max;
dPhi_des  = amp*(3*tanhT2.^4 - 4*tanhT2.^2 + 1)/Phi_0_max;
ddPhi_des = amp*(-4*tanhT2.*(3*tanhT2.^4 - 5*tanhT2.^2 + 2))/Phi_0_max;

% the following "if 1" or "if 0" construction is useful when you'd
%   like to turn things on and off at different times

if 1
 figure
  plot(T0, Phi_des*180/pi)
  grid on, zoom on
 title('\phi')
end

% if the Phi_des curve satisfies  | phi(t) | < pi/2, all t,
%   then it can be part of a feasible trajectory
%   and we can determine the control portion of that traj

% params given in dynamics.c
GG = 9.81;  LL = 0.5;

% from dynamics.c
%    dphid = G_L*s_phi - o_L*c_phi*uu;
% rewrite as
%    G s_phi - L ddphi = c_phi u

if  max( abs(Phi_des) ) < 60*pi/180   % be somewhat conservative!
  Udes = ( GG*sin(Phi_des) - LL*ddPhi_des ) ./ cos(Phi_des);
else
  Udes = 0*T0;  % or whatever you might *invent* (??)
end

Xdes = [ Phi_des dPhi_des ];

% note: this is but one of *many* ways to
%   get *started* with playing with system trajectories
% DO NOT get stuck thinking that it is the *right* way !!!





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




% let's now try using a constant (time-invariant) state feedback
%   to attempt tracking some trajectories and some non-trajectories

% design a TI LQR feedback K about
%   the trivial traj (inverted eq pt)  (x,u)(t) == ( (0,0), 0 )

% get A & B

% look at dynamics_m.c to see the calling sequence
% [dx,y, fxu_x,fxu_u, q_fxu_x_x,q_fxu_x_u,q_fxu_u_u] = dynamics_m(x, u, wt, q);
[~,~,A,B] = dynamics_m([0 0], 0);

% choose the naive choice Q = I, R = 1
[K,~,E] = lqr(A,B,diag([1 1]),1);

% it is often helpful to look at the 'initial condition response'
%   of the linear approximation
% to see whether the constant K looks good to you
% this can be done using a SS (state space) system

sys = ss(A-B*K, B, [ [1 0]; -K ], 0 );
% here, I've specified outputs as  phi and u
%   (and cheated on the D matrix, it's 2x1 ...)

if 0
 % play with a few responses
 figure, initial(sys, [1 0]), grid on, zoom on
 figure, initial(sys, [0 1]), grid on, zoom on
 figure, initial(sys, [1 1]), grid on, zoom on
end
 


% PLAY a bit HERE!  (with different K's and different IC's)
% now let's get one of the linear responses
%   and compare it with the *nonlinear* response

% set up nonlin stuff

% get system sizes
[NS, NI, NO, NW, NWL, NWC] = sys_sizes_m;

NS2 = NS*NS; NN12 = NS*(NS+1)/2;


% Xi, Ui  will be  the trajectory to (try to) track [Xdes,Udes is for cost]
%   here, read Xi as  x sub i (x_i)  [ rather than  \xi ]
%   the 'sub i' derives from it's use in the optimization iteration
Xi = zeros( length(T0), NS );
Ui = zeros( length(T0), NI );

Wt = [ ];           % *no* exogenous input to dynamics
Wlt = [Xdes Udes];  % exogenous input to incremental cost



% PLAY with LQR design !!
Qr = diag([1 1]);
Rr = 1;

[K,~,E] = lqr(A, B, Qr, Rr);

Kr = repmat(K,length(T0),1);    % Kr(t) for nonl sim



% PLAY with initial condition !!
x00 = [1 0];
x00 = [0.1 0];

%%Exercise%%
% determine (for some choices of Q,R),
%   the angle phi0 (to 0.01 degree)
% such that
%   the nonlinear response from
%     x00 = [ phi0 0 ]
%   converges to [0 0]
% but the response from
%     x00 = [ phi0 0] + [ 0.01*pi/180 0]
%   does not.
% ( [phi0 0] is clsoe to the boundary of the region of attraction )
% [ making the region of attraction large may or may *not* be desireable for us ]

phi0 = 66.7*pi/180;   x00 = [phi0 0];


% nonl_K implements the projection operator (+ cost evaluation)
%   note the (extended) input

ip_opts = simset(opts, 'InitialState', [x00 0]);
[Ti, aXi, aYi] = sim('nonl_KS', T0, ip_opts, [T0 [Xi Ui] Kr Wt Wlt]);
     %  the 'a' in aXi, aYi means 'augmented' (why?)

Xnonl = aXi(:,1:NS);
Unonl = aYi(:,1:NI);

% use the same time grid for *both* lin and nonlin
[Ylin,Tlin,Xlin] = initial(sys, x00, T0);
Ulin = Ylin(:,2);

% now, compare
%   I've found them close (using Qr=I,Rr=1)
%     when  x00 = [0.1 0]
%   and not so close when  x00 = [1 0]

figure
  plot(T0, [Xnonl(:,1), Xlin(:,1)]*180/pi)
  grid on, zoom on
title('\phi(t), nonlinear and linear responses')

if 0    % I don't *always* care to see this ...
 figure
  plot(T0, [Unonl, Ulin])
  grid on, zoom on
 title('u(t), nonlinear and linear responses')
end



%%Exercise%%
% try to get a feel for how the
%   convergence to eq pts (various eq pts)
% behavior changes
%   when you choose different constant Q,R
%
% sometimes, one can find one choice such that
%   the transient behavior looks similar (and acceptable)
%   over an intereting region of the equilibrium manifold
%
% such a choice can lead to good time-varying regulator performance
%   for the projection operator



% I sometimes use a 'return' to allow me to partially execute a matlab script
%  e.g., to get started "all at once"

% return

%************************


% here, I include calls for the TV regulator design 

% less play comments ...



x00 = Xdes(1,:);
u00 = Udes(1,:);

% use LQR to design a P1 for terminal cost
%   ... hopefully, we can also use TI K(.) for simple tracking

% Q and R defined by cost_params_h.m

% linearize about terminal state
x11 = Xdes(end,:);
u11 = Udes(end,:);
[~,~,A1,B1] = dynamics_m(x11,u11,[]);

% Qreg = diag([1 1]);  % change this as needed
% Rreg = 1;
%
% set Qreg and Rreg in QR_params.h  (and mex)
% and grab those
[Qreg, Rreg] = QR_params_m;

[Kf,Pf,Ef] = lqr(A1,B1,Qreg,Rreg);

% organize the lower triangle of Pf in a flat manner
Pf_ = [];
for ind=1:NS
  Pf_ = [Pf_ Pf(ind,1:ind)];
end

% pretend that Xdes,Udes is a trajectory
Xi = Xdes;
Ui = Udes;

% OR, use a trivial trajectory
% Xi = repmat(x00,length(T0),1);
% Ui = repmat(u00,length(T0),1);


% integrate Riccati eqn *backward* in time --- note flipud

Popts = simset(opts, 'InitialState', Pf_);
[Tb,Pb,Kb] = sim('LQR_KrS',T0,Popts, [T0 flipud([Xi Ui])]);

Kr = flipud(Kb(:,1:(NI*NS)));
Pi = flipud(Pb(:,1:NN12));

if 0
 figure
  plot(T0,[Pi repmat(Pf_,length(T0),1)])
  grid on, zoom on
 title('P')

 figure
  plot(T0,[Kr repmat(Kf,length(T0),1)])
  grid on, zoom on
 title('K')
end

% now use nonl_K to try tracking Xi,Ui
Xi = Xdes;
Ui = 0*Udes;

ip_opts = simset(opts, 'InitialState', [x00 0]);
[Ti, aXi, aYi] = sim('nonl_KS', T0, ip_opts, [T0 [Xi Ui] Kr Xdes Udes]);


if 1
 figure
  plot(T0,[Xdes(:,1) aXi(:,1)]*180/pi)
  grid on, zoom on
 title('\phi')

 figure
  plot(T0,[Udes aYi(:,1)]*180/pi)
  grid on, zoom on
 title('u')
end


