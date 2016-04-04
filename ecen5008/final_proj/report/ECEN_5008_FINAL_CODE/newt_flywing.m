% newt_slcar.m
%
%   for the flying wing system,
%     find system trajectory minimizing L_2 trajectory error
%
%  John Hauser
%  nov/dec 06
%  sept 09 stuttgart
%  oct 15 boulder
%  nov 15 boulder



% verbose provides the capability
%   to control some printing, plotting, and pausing options
%
%   verbose = 0;  % (pretty) quiet run, no plots until end
%   verbose = 1;  % more text, no plots until end
%   verbose = 2;  % plots (x,u)_i, (z,v)_i at each iteration
%   verbose = 3;  % also plots line search (extra evaluations)
%   verbose = 4;  % plot everything and pause after search direction

% verbose = 4;
verbose = 3;
% verbose = 2;
% verbose = 1;
% verbose = 0;


% get system sizes: NS, NI, NO, etc.
[NS, NI, NO, NW, NWL, NWC] = sys_sizes_m;

NS2 = NS*NS;
NN12 = NS*(NS+1)/2;


% options for 'sim'
opts = simset('Solver','ode45', ...
	      'AbsTol',1e-8,'RelTol',1e-6); %, ...


% turn off some warnings that are not applicable to our MDLs
%   also, throw an error if Riccati wants to blow up in Prq2_KvS
load_system('LQR_KrS');
load_system('nonl_KS');
load_system('Prq1_KvS');
load_system('Prq2_KvS');
load_system('linK_zvS');
set_param('LQR_KrS', 'InitInArrayFormatMsg','None');
set_param('nonl_KS', 'InitInArrayFormatMsg','None');
set_param('Prq1_KvS','InitInArrayFormatMsg','None');
set_param('Prq2_KvS','InitInArrayFormatMsg','None','MinStepSizeMsg','error');
set_param('linK_zvS','InitInArrayFormatMsg','None');
save_system('LQR_KrS');
save_system('nonl_KS');
save_system('Prq1_KvS');
save_system('Prq2_KvS');
save_system('linK_zvS');


fig_base = 1000;


% get screen size for plotting
scrsz = get(0,'ScreenSize');
h = figure;
pos  = get(h,'Position');
opos = get(h,'OuterPosition');
close(h);

wscr = scrsz(3);
hscr = scrsz(4);

wfig = pos(3);
hfig = pos(4);
hdfig = opos(4) -  pos(4);
tpfig = opos(2) + opos(4);

wnfig = round( (wscr - wfig)/2 - 5 );
hnfig = round( hfig/wfig*wnfig );
ytnfig = tpfig - hdfig - hnfig;
ybnfig = ytnfig - hdfig - hnfig - 5;

XiPos = [0 ytnfig wnfig hnfig];
UiPos = [0 ybnfig wnfig hnfig];
ZiPos = [wscr-wnfig ytnfig wnfig hnfig];
ViPos = [wscr-wnfig ybnfig wnfig hnfig];
DescPos = [pos(1)+wfig-wnfig ybnfig wnfig hnfig-(hfig-hnfig)];




% save des_traj  Tdes Xdes Udes  T0 X0 U0
load des_traj  % get Tdes Xdes Udes  T0 X0 U0
T0 = Tdes;

% for use in dynamics and cost
Wt  = [ ];  % no time dependent input to dynamics
Wlt = [ Xdes Udes ];  % incremental cost for L2 traj error


% let's suppose beginning and end are (close to being) equilibrium points
x00 = X0(1,:);
u00 = U0(1,:);

x11 = Xdes(end,:);
u11 = Udes(end,:);


% to get started, let's design a K around straight running
% x0 = [ x00(1,1) 0 0 ];
% u0 = [ 0 0];
x0 = Xdes(1,:);
u0 = Udes(1,:);

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

% for now, try a time invariant K # Initial K for optimization?
[K,~,E] = lqr(A0, B0, Qr, Rr);

% since NI > 1, we need to flatten K
K_ = [];
for ii=1:size(K,1)
  K_ = [ K_, K(ii,:) ];
end

Kr = repmat(K_,length(T0),1);    % Kr(t) for nonl sim

if 0
 % we suppose that x00,u00 is (close to being) an equilibirum point
 X0 = repmat(x00,length(T0),1);
 U0 = repmat(u00,length(T0),1);

 Xi = X0;
 Ui = U0;
end

if 0
 % pretend that Xdes,Udes *is* a trajectory ...
 X0 = Xdes;
 U0 = Udes;

 Xi = X0;
 Ui = U0;
end

Xi = X0;
Ui = U0;

% BUT, we'll try to get an actual trajectory
nonl_opts = simset(opts, 'InitialState', [x00 0]);
[Ti, aXi, aYi] = sim('nonl_KS', T0, nonl_opts, [T0 [Xi Ui] Kr Wt Wlt]);

Xi = aXi(:,1:NS);
Ui = aYi(:,1:NI);

X0 = Xi;
U0 = Ui;


% get linearization about the final desired condition
%   for terminal Riccati P(T) cost & regulator  [heuristics]
[~,~,A1,B1] = dynamics_m(x11,u11);


% get params for regulator and cost

% for terminal cost
[~, ~, ~, Q, ~, R] = cost_m(x11, u11, [ x11 u11 ]);
% [lxu, lxu_x, lxu_u, lxu_x_x, lxu_x_u, lxu_u_u] = cost_m(x, u, wlt);

% get terminal cost P1  (call it Pf and Pf_)  [heuristic]
[Kf,Pf,Ef] = lqr(A1,B1,Q,R);


% for LQR
[Qreg, Rreg] = QR_params_m;

% get regulator terminal P  [heuristic]
%    (just fix it for desired terminal equil state)
[Kreg,Preg,Ereg] = lqr(A1,B1,Qreg,Rreg);



% organize the lower triangle of Pf (and Preg) in a flat manner
Pf_ = [];
Preg_ = [];
for ind=1:NS
  Pf_ = [Pf_ Pf(ind,1:ind)];
  Preg_ = [Preg_ Preg(ind,1:ind)];
end


% always start nonl_K at the *same* init cond
nonl_opts = simset(opts, 'InitialState', [x00 0]);

% regulator also always starts at same init (final) cond
Kr_opts = simset(opts, 'InitialState', Preg_);


% in case we are silly and try an open loop initial trajectory
if 0
  Kr = zeros(length(T0),NS*NI);

  [Ti, aXi, aYi] = ...
     sim('nonl_KS', T0, ip_opts, [T0 [Xi Ui] Kr Wt Wlt]);

  X0 = aXi(:,1:NS);
  U0 = aYi(:,1:NI);
end

if 0
  figure,plot(T0,X0),grid on, zoom on
  pause
end


% x sub i, etc
Xi = X0;
Ui = U0;


% algorithm params
gam_alf = 0.4;	gam_beta = 0.7;
% MaxIter = 40;	DescentTol = 1e-6;
MaxIter = 40;	DescentTol = 1e-7;
max_gam = 1;

Colors = 'bgrcmyk';	nColors = 7;	% for plotting

% to store results

Descent = [];
Quad = [];
Vals = [];
MinCost = [];
MaxZV = [];
Method = [];
Gammas = [];
iGammas = [];
CompTimes = [];


for ind = 0:MaxIter

  compTime = 0;

  %************************************************
  % evaluate the functional (after designing K)
  %************************************************

  % design Kr for Xi,Ui

  cpu = cputime;

  % integrate Riccati eqn *backward* in time --- note flipud
  Tb = T0(end) - flipud(T0);
  [Tb,Pb,Kb] = sim('LQR_KrS',Tb,Kr_opts, [Tb flipud([Xi Ui])]);
  Kr = flipud(Kb(:,1:(NI*NS)));

  % project and evaluate

  [Ti, aXi, aYi] = ...
    sim('nonl_KS', T0, nonl_opts, [T0 [Xi Ui] Kr Wt Wlt]);

  if 0
    max_dXU = max(abs( [ aXi(:,1:NS)-Xi, aYi(:,1:NI)-Ui ] ));
    fprintf('max_dXU:');
    fprintf(' %g', max_dXU);
    fprintf('\n');
  end

  Xi = aXi(:,1:NS);
  Ui = aYi(:,1:NI);
  if NO > 0
    Yi = aYi(:,NI+1:NI+NO);
  end

  Xf = Xi(end,:);
  dXf = Xf - x11;
  val = aXi(end,NS+1) + 1/2*dXf*Pf*dXf';
  %  Ji = val - aXi(:,NS+1);   % cost to go

  compTime = compTime + cputime - cpu;


  if verbose > 1
    figure(fig_base + 1), set(gcf,'Position',XiPos);
      plot(T0,Xdes,'-.',T0,Xi,'-')
      grid on, zoom on
    title(sprintf('X%d',ind))
 
    figure(fig_base + 2), set(gcf,'Position',UiPos);
      plot(T0,Udes,'-.',T0,Ui,'-')
      grid on, zoom on
    title(sprintf('U%d',ind))
  end


% pause

  %************************************************
  % get descent direction  zeta = (z,v) 
  %   and compute  Dg(xi).zeta          (= descent)
  %         and    D2g(xi).(zeta,zeta)  (= quad)
  %************************************************

  Prq0 = [Pf_ dXf*Pf dXf*Pf];
  Popts = simset(opts, 'InitialState', Prq0);

  cpu = cputime;

  % heuristic for using FIRST or SECOND order method
  % if ind<3 | descent_toggle < 1,
  if ind<2      %   | descent_toggle < 1,

    Tb = T0(end) - flipud(T0);
    [Tb,Prq,Kvb] = ...
        sim('Prq1_KvS',Tb,Popts, ...
            [Tb flipud([Xi Ui Kr Wt Wlt])]);
    Method=[Method;1];

  else

    % try using second order method
    try
      Prq_cpu = cputime;
      Tb = T0(end) - flipud(T0);
      [Tb,Prq,Kvb]=sim('Prq2_KvS',Tb,Popts, ...
          [Tb flipud([Xi Ui Kr Wt Wlt])]);
      Prq_cpu = cputime - Prq_cpu;
      % fprintf('Prq2_cpu = %g\n',Prq_cpu);

      Method = [Method;2];


    catch
      % disp(lasterr);
      % disp('switching to first order descent');
      Prq_cpu = cputime;
      Tb = T0(end) - flipud(T0);
      [Tb,Prq,Kvb]=sim('Prq1_KvS',Tb,Popts, ...
          [Tb flipud([Xi Ui Kr Wt Wlt])]);
      Prq_cpu = cputime - Prq_cpu;
      % fprintf('Prq1_cpu = %g\n',Prq_cpu);

      Method=[Method;1];

    end

  end

  Ki = flipud(Kvb(:,1:(NI*NS)));
  Vo = flipud(Kvb(:,((NI*NS)+1):((NI*NS)+NI)));
  Ri = flipud(Prq(:,(NN12+1):(NN12+NS)));
  Qi = flipud(Prq(:,(NN12+NS+1):(NN12+NS+NS)));
  Pi = flipud(Prq(:,1:NN12));

  [Tb,aZi,aVi]=sim('linK_zvS',T0,opts, ...
                  [T0 Xi Ui Ki Vo Qi Wt Wlt]);
  Zi = aZi(:,1:NS);
  Vi = aVi(:,1:NI);

  Zf = Zi(end,:);
  descent = aZi(end,NS+1) + dXf*Pf*Zf';
  quad = aZi(end,NS+2) + Zf*Pf*Zf';

  compTime = compTime + cputime - cpu;

  if verbose > 0
    fprintf('ind: %d: %g [%d]\n', ind, val,Method(end));
    if verbose > 1
      fprintf('  descent: %g   term: %g\n', descent, dXf*Pf*Zf');
      fprintf('  quad: %g   term: %g\n', quad, 1/2*Zf*Pf*Zf');
    
    end
  end

  Vals = [Vals; val];
  Descent = [Descent; descent];
  Quad = [Quad; quad];

  % restrict the step length (max_gam)
  %   to some max change in phi, theta (say, 30 degrees) ***special***
  maxZi = max(max(abs(Zi(:,[1]))));

  % MinCost = [MinCost; min_cost];
  MaxZV = [MaxZV; [ maxZi, max(abs(Vi))] ];

  if verbose > 1
    figure(fig_base + 3), set(gcf,'Position',ZiPos);
      plot(Ti, Zi)
      grid on, zoom on
     title(sprintf('Z%d',ind))

    figure(fig_base + 4), set(gcf,'Position',ViPos);
      plot(Ti, Vi)
      grid on, zoom on
    title(sprintf('V%d',ind))

    if ind>0
     figure(fig_base + 5), set(gcf,'Position',DescPos);
      plot(0:(length(Descent)-1), log10(-Descent),'g-', ...
           0:(length(Descent)-1), log10(-Descent),'r+', ...
           0:(length(Vals)-2), log10(-diff(Vals)),'c-o')
      grid on, zoom on
     title('log10 -Descent')
    end

  end

  if verbose > 3
    pause
  end

  if ind > 0
    % save trace of computation
    eval(sprintf('T%i = T0;',ind));
    eval(sprintf('X%i = Xi;',ind));
    eval(sprintf('U%i = Ui;',ind));
  end

  if descent > -DescentTol,
    fprintf('convergence achieved: descent: %g [%d]\n', descent,Method(end));
    break;
  end

  %************************************************
  % backtracking (armijo) line search
  %************************************************
  %
  %   using  gam_alf in ( 0, 1/2 ),  gam_beta in ( 0, 1 )
  % 
  %   accept gamma = gam_beta ^ k, k = 0, 1, ...
  %
  %     g(xi + gamma*zeta) < g(xi) +  gam_alf * descent * gamma
  %

  gam_accept = 0;

  % ***special***
  %max_gam = 0.5/maxZi;     % gam=1 corresponds to 0.5 rads in phi, theta max
  %max_gam = 1.0/maxZi;     % gam=1 corresponds to 1.0 rads in phi, theta max
  max_gam = 1.0; %***************

  gammai = min(1.0,max_gam);
  % gammai = 1.0;
  gamsi = [];
  gCost = [];

  cpu = cputime;

  for k = 0:12,
    [Ti,aXi,aYi] = ...
      sim('nonl_KS',T0,nonl_opts, ...
          [T0 ([Xi Ui]+gammai*[Zi Vi]) Kr Wt Wlt]);
    g_dXf = aXi(end,1:NS) - x11;
    gamsi = [gamsi; gammai];
    gCost = [gCost; aXi(end,NS+1) + 1/2*g_dXf*Pf*g_dXf'];

    if verbose > 1
      fprintf(' gCost: %g\n', gCost(end));
    end

    if gCost(end) < val  +  gam_alf * descent * gammai,
      gam_accept = 1;
      % igam = find(gCost == min(gCost));
      igam = find(gCost == gCost(end));  % silly quick fix
      iGammas = [iGammas; igam];
      gammai = gamsi(igam);
      Gammas = [Gammas; gammai];
      break;
    end

    gammai = gam_beta*gammai;
  end

  compTime = compTime + cputime - cpu;

  if verbose > 1
    gammaS = -descent/quad;
    fprintf('  gammaS:  %g   predict:  %g\n', gammaS, gCost(1)-descent^2/quad/2);
  end

  %*******************************

  if verbose > 2
    gams = (0:.1:1.2) * min(1.0,max_gam);
    Cost = zeros(size(gams));
    len2 = ceil(3/4*length(gams));

    for k = 1:length(gams)
      gamma = gams(k);
      [Ti,aXi,aYi] = ...
        sim('nonl_KS',T0,nonl_opts, ...
            [T0 ([Xi Ui]+gamma*[Zi Vi]) Kr Wt Wlt]);
      g_dXf = aXi(end,1:NS) - x11;
      Cost(k) = aXi(end,NS+1) + 1/2*g_dXf*Pf*g_dXf';
    end
  
    figure(fig_base + 10 + ind)
      plot(gams,min(Cost,2*Cost(1)), ...
  	gamsi, min(gCost, 2*Cost(1)), 'x', ...
  	... % [0 .5], Cost(1)+descent*[0 .5], ...
    	... % [0 1.3], Cost(1)+descent/2*[0 1.3], ...
  	gams(1:len2), Cost(1)+descent*gams(1:len2), ...
    	gams, Cost(1)+descent/2*gams, ...
    	gams, Cost(1) + descent*gams + quad*gams.*gams/2, ...
    	gams, Cost(1)+0.4*descent*gams ...
	... % [0 1.3], Cost(1)+0.4*descent*[0 1.3] ...
          )
      title(sprintf('order %d: [X%d U%d] + gamma*[Z%d V%d]',Method(end),ind,ind,ind,ind));
      grid on, zoom on
    drawnow
  end

  if ( ~ gam_accept ),
    % error( ) ??
    disp('*** unable to satisfy step size --- aborting');
    break;
  end

  if verbose > 1
    fprintf('  gammai = %g\n', gammai);
  end


  %************************************************
  % update
  %************************************************

  cpu = cputime;

  [Ti,aXi,aYi] = ...
    sim('nonl_KS',T0,nonl_opts, ...
        [T0 ([Xi Ui]+gammai*[Zi Vi]) Kr Wt Wlt]);

  Xf = aXi(end,1:NS);
  g_dXf = aXi(end,1:NS) - x11;
  val = aXi(end,NS+1) + 1/2*g_dXf*Pf*g_dXf';

  Xi = aXi(:,1:NS);
  Ui = aYi(:,1:NI);
  Ji = val - aXi(:,NS+1);

  compTime = compTime + cputime - cpu;

  if verbose > 1
    fprintf('  val: %g   term: %g   norm(dXf): %g\n', val, 1/2*g_dXf*Pf*g_dXf', norm(g_dXf));
    fprintf('  Xf: %g %g %g\n', Xf(1:NS));
  end

  % Xf = Xi(end,:);
  % dXf = Xf - x11;

  if 0    % verbose>1,
    figure(fig_base + 1)
      plot(To,Xo(:,[1]),'-.',T0,X0(:,[1]),'-.',Ti,Xi(:,[1]),'-')
      grid on, zoom on
    title('Xo and Xi')
 
    figure(fig_base + 2)
      plot(To,Uo,'-.',T0,U0,'-.',Ti,Ui,'-')
      grid on, zoom on
    title('Uo and Ui')
  end

  CompTimes = [CompTimes; compTime];

end

fprintf('newt CompTime: %g\n', sum(CompTimes));


Xopt = Xi;
Uopt = Ui;
Kopt = Ki;


fig_base = fig_base + 1000;

figure(fig_base + 1), set(gcf, 'Position', XiPos+[10 -10 0 0]);
  plot(T0,Xdes(:,1),'m-.')
hold on
  for k=0:(length(Vals)-1),
    eval(sprintf('plot(T%d,X%d(:,1),''%c'');',k,k,Colors(1+mod(k,nColors))))
  end
hold off
  grid on, zoom on
title('X1i vs T')

figure(fig_base + 2), set(gcf, 'Position', XiPos+2*[10 -10 0 0]);
  plot(T0,Xdes(:,2),'m-.')
hold on
  for k=0:(length(Vals)-1),
    eval(sprintf('plot(T%d,X%d(:,2),''%c'');',k,k,Colors(1+mod(k,nColors))))
  end
hold off
  grid on, zoom on
title('X2i vs T')

figure(fig_base + 3), set(gcf, 'Position', XiPos+3*[10 -10 0 0]);
  plot(T0,Xdes(:,3),'m-.')
hold on
  for k=0:(length(Vals)-1),
    eval(sprintf('plot(T%d,X%d(:,3),''%c'');',k,k,Colors(1+mod(k,nColors))))
  end
hold off
  grid on, zoom on
title('X3i vs T')

figure(fig_base + 4), set(gcf, 'Position', XiPos+4*[10 -10 0 0]);
  plot(T0,Xdes(:,4),'m-.')
hold on
  for k=0:(length(Vals)-1),
    eval(sprintf('plot(T%d,X%d(:,4),''%c'');',k,k,Colors(1+mod(k,nColors))))
  end
hold off
  grid on, zoom on
title('X4i vs T')


figure(fig_base + NS+1), set(gcf, 'Position', UiPos+1*[10 -10 0 0]);
  plot(T0,Udes(:,1),'m-.')
hold on
  for k=0:(length(Vals)-1),
    eval(sprintf('plot(T%d,U%d(:,1),''%c'');',k,k,Colors(1+mod(k,nColors))))
  end
hold off
  grid on, zoom on
title('U1i vs T')

figure(fig_base + NS+2), set(gcf, 'Position', UiPos+2*[10 -10 0 0]);
  plot(T0,Udes(:,2),'m-.')
hold on
  for k=0:(length(Vals)-1),
    eval(sprintf('plot(T%d,U%d(:,2),''%c'');',k,k,Colors(1+mod(k,nColors))))
  end
hold off
  grid on, zoom on
title('U2i vs T')

figure(fig_base + NS + NI + 1), set(gcf,'Position',DescPos);
  plot(0:(length(Descent)-1), log10(-Descent),'g-', ...
       0:(length(Descent)-1), log10(-Descent),'r+', ...
       0:(length(Vals)-2), log10(-diff(Vals)),'c-o')
  grid on, zoom on
title('log10 -Descent')

v = @(t) interp1(T,X1(:,1),t);
gamma = @(t) interp1(T,X1(:,2),t);


%%%% Integrate to get x and z
vdes = @(t)interp1(T,Xdes(:,1),t);
gammades = @(t) interp1(T,Xdes(:,2),t);
kinemat = @(t,x) [v(t)*cos(gamma(t)); -v(t)*sin(gamma(t))];
kinematdes = @(t,x) [vdes(t)*cos(gammades(t)); -vdes(t)*sin(gammades(t))];
[Ti,XZ] = ode45(kinemat,T,[0 0],odeset('AbsTol',1e-8,'RelTol',1e-6));
[Tides, XZdes] = ode45(kinematdes,T,[0 0],odeset('AbsTol',1e-8,'RelTol',1e-6));
figure,plot(XZdes(:,1),-XZdes(:,2),XZ(:,1),-XZ(:,2)),grid on,zoom on,axis equal,title('Optimized Wing Path')
xlabel('X  (m)');
ylabel('-Z (m)');
legend('Desired', 'Actual')

% return
