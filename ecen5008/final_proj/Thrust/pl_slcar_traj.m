% plot slcar trajectory
%
%   T0 Xi Ui Yi Xdes Udes Alat_des   should exist
%
% JH
% nov 15  boulder

figure, plot(T0, Xdes(:,1), '-.', T0, Xi(:,1)), grid on, zoom on, title('V')
figure,plot(T0,Xdes(:,2)*180/pi,'-.',T0,X1(:,2)*180/pi),grid on,zoom on,title('\beta')
figure,plot(T0,Xdes(:,3)*180/pi,'-.',T0,X1(:,3)*180/pi),grid on,zoom on,title('\omega')

figure,plot(T0,Alat_des,'-.',T0,Yi(:,2)),grid on,zoom on,title('a_l_a_t')

figure,plot(T0,Udes(:,1),'-.',T0,Ui(:,1)),grid on,zoom on,title('u_1')
figure,plot(T0,Udes(:,2),'-.',T0,Ui(:,2)),grid on,zoom on,title('u_2')


% use anonymous functions with ode45 to integrate the kinematics
v = @(t) interp1(T0,Xi(:,1),t);
beta = @(t) interp1(T0,Xi(:,2),t);
omega = @(t) interp1(T0,Xi(:,3),t);
kinemat = @(t,x) [ v(t)*cos(x(3)+beta(t)); v(t)*sin(x(3)+beta(t)); omega(t) ];

psi0 = -65*pi/180;
[Ti,XYPsi] = ode45(kinemat,T0,[0 0 psi0],odeset('AbsTol',1e-8,'RelTol',1e-6));

figure,plot(XYPsi(:,2),XYPsi(:,1)),grid on,zoom on,axis equal,title('ground track')
