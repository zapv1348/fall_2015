% This .m-file circular_current computes ray paths through a current field. As 
% an example, we use the circular current example in the modeling book.
close all;

% Define the computational domain
n = 80;                                % Number of grid intervals in x- and y-
                                       % directions used to represent the current
de = 200e3;                            % field over a physical domain extending
[Y,X] = ndgrid(linspace(-de,de,n+1));  % X and Y coordinates at node points

% Set up the circular current field U, V
ri = 40e3; ro = 160e3;                 % Edges of current at dist. ri and ro resp.
rc = (ri+ro)/2; 
sp = 2;                                % Current speed 2 m/s.
U = zeros(n+1); V = zeros(n+1);        % Reserve space for the velocity fields
av = zeros(n+1);
r       = sqrt(X.^2+Y.^2);             % Find distance to center and note 
cur     = r>ri & r<ro;                 % where the current is present

% Plot the circular current
ang = linspace(0,2*pi,101); si = sin(ang)*1e-3; co = cos(ang)*1e-3;
plot(ri*si,ri*co,':',rc*si,rc*co,':',ro*si,ro*co,':')
axis ([-200 200 -200 200]); axis square; hold on;

% Calculate velocities in the current (and elsewhere)
av(cur) = sp*4*(r(cur)-ri).*(r(cur)-ro)./(r(cur)*(ro-ri)^2);
U(cur)  = -av(cur).*Y(cur);            % Velocity in the x-direction
V(cur)  =  av(cur).*X(cur);            % Velocity in the x-direction
h = 2*de/n;                            % Grid spacing
[UX,UY] = gradient(U,h); [VX,VY] = gradient(V,h);   % Calculate UX UY VX VY 

% Prepare for 2-D interpolation during the time stepping
Iu  = griddedInterpolant(Y,X,U, 'cubic');   
Iv  = griddedInterpolant(Y,X,V, 'cubic');   
Iux = griddedInterpolant(Y,X,UX,'cubic');
Iuy = griddedInterpolant(Y,X,UY,'cubic');
Ivx = griddedInterpolant(Y,X,VX,'cubic');
Ivy = griddedInterpolant(Y,X,VY,'cubic');

% Define physical parameters; needed for initial k1 and k2 values 
g = 9.8;                               % Acceleration of gravity
t = 10;                                % Time period of incoming waves (in seconds)
s = 2*pi/t;                      	   % Angular time frequency
k = s^2/g;                             % Magnitude of initial wave vector

% Loop over the rays
xykk = zeros(4,1);  
dt = 1000;                             % Specify time step (in seconds).
for lt = -360:10:360
   xykk = [max( lt-180,-180)*1e3; ...  % Give start x- and y-locations for the   
           max(-lt-180,-180)*1e3; ...  % rays shown in Figure I2.4-3.
           k/sqrt(2); ...              % Component in x-dir (45 degree)
           k/sqrt(2)];                 % Component in y-dir (45 degree)
   
   % Run ODE for ray until it exits the domain, using 2-stage second order RK
   ct = 1; clear x y;
   while 1                             % Keep stepping until stop criterion      
      rhs_xykk = c_rhs(xykk,              X,Y,Iu,Iv,Iux,Iuy,Ivx,Ivy); % RK stage 1 ; evaluate RHS
      rhs_xykk = c_rhs(xykk+dt/2*rhs_xykk,X,Y,Iu,Iv,Iux,Iuy,Ivx,Ivy); % RK stage 2 ; evaluate RHS
      if max(xykk(1:2))>X(1,end-2)||min(xykk(1:2))<X(1,1); break;end % Stop when ray exits domain
      xykk  = xykk + dt*rhs_xykk;      % Update vector with x,y,k1,k2            
      x(ct) = xykk(1); y(ct) = xykk(2);% Store x,y-values to plot the ray.       
      ct = ct+1;
   end
                                       % Plot the ray;
   plot(x(:)*1e-3,y(:)*1e-3,'k-')      % Scale down by 1e-3 to get unit km.   
end
