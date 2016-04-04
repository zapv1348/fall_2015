
ta = 5;
tb = 7.5;
tc = 12.5;
td = 15;
tf = 20;

T0 = 0:0.01:ta;
T1 = (ta+0.01):0.01:tb;
T2 = (tb+0.01):0.01:tc;
T3 = (tc+0.01):0.01:td;
T4 = (td+0.01):0.01:tf;

gamma1 = pi*ones(1,length(T0));
gamma2 = pi+((pi/4)/(tb-ta))*(T1-ta);
gamma3 = (5*pi/4)*ones(1,length(T2));
gamma4 = (5*pi/4)-((pi/4)/(td-tc))*(T3-tc);
gamma5 = pi*ones(1,length(T4));

T = [T0 T1 T2 T3 T4];
Gamma = [gamma1 gamma2 gamma3 gamma4 gamma5];
plot(T,Gamma)

