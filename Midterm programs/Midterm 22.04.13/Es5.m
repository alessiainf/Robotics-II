clc; close all; clear;
syms q1 q2 a b c l1 l2 real

q=[q1 q2];
qsub=[1 pi/2];

M=[a, b*cos(q2);
    b*cos(q2), c];
M=subs(M, q, qsub)

p=[l2*cos(q2);
    q1+l2*sin(q2)];
J=jacobian(p, q);
J=subs(J, q, qsub)

%taua=....=[a*yd_dd+g1(q);
%               0];

%taub=....=taua;