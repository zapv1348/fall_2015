function [ uvnew ] = RK4( u1, u2, u3, v1, v2, h, Pr )
% uses the functions F1-F3, G1, G2 to perform one iteration of RK4 with
% step size h, returns a vector with updated values for u1, u2, u3, v1, v2
k11 = h * F1(u1, u2, u3);
k12 = h * F2(u1, u2, u3);
k13 = h * F3(u1, u2, u3);
k21 = h * F1(u1 + k11/2, u2 + k12/2, u3 + k13/2);
k22 = h * F2(u1 + k11/2, u2 + k12/2, u3 + k13/2);
k23 = h * F3(u1 + k11/2, u2 + k12/2, u3 + k13/2);
k31 = h * F1(u1 + k21/2, u2 + k22/2, u3 + k23/2);
k32 = h * F2(u1 + k21/2, u2 + k22/2, u3 + k23/2);
k33 = h * F3(u1 + k21/2, u2 + k22/2, u3 + k23/2);
k41 = h * F1(u1 + k31, u2 + k32, u3 + k33);
k42 = h * F2(u1 + k31, u2 + k32, u3 + k33);
k43 = h * F3(u1 + k31, u2 + k32, u3 + k33);
u1new = u1 + (k11 + 2*k21 + 2*k31 + k41)/6;
u2new = u2 + (k12 + 2*k22 + 2*k32 + k42)/6;
u3new = u3 + (k13 + 2*k23 + 2*k33 + k43)/6;
l11 = h * G1(u1, v1, v2, Pr);
l12 = h * G2(u1, v1, v2, Pr);
l21 = h * G1(u1, v1 + l11/2, v2 + l12/2, Pr);
l22 = h * G2(u1, v1 + l11/2, v2 + l12/2, Pr);
l31 = h * G1(u1, v1 + l21/2, v2 + l22/2, Pr);
l32 = h * G2(u1, v1 + l21/2, v2 + l22/2, Pr);
l41 = h * G1(u1, v1 + l31, v2 + l32, Pr);
l42 = h * G2(u1, v1 + l31, v2 + l32, Pr);
v1new = v1 + (l11 + 2*l21 + 2*l31 + l41)/6;
v2new = v2 + (l12 + 2*l22 + 2*l32 + l42)/6;
uvnew = [u1new, u2new, u3new, v1new, v2new];


end

