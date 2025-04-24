clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real

%the trajectory is:
syms a b k t T
qd=[a+b*(1-cos((pi*t)/T)); k];
qd_dot=diff(qd, t);

%inertia matrix
M = [m1*dc1^2 + m2*l1^2 + m2*q2^2 + I1 + I2, -l1*m2;
                                -l1*m2,     m2];

M=subs(M, [q1,q2], [qd(1), qd(2)]);

p=M*qd_dot;

h = int(p, t, 0, T)