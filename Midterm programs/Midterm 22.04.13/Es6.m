clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 real
syms q1_dot q2_dot q3_dot real
syms a b c real
syms rc1x rc1y rc2x rc2y rc3x rc3y real

M=[a, b*cos(q2);
    b*cos(q2), c];

q=[q1, q2];

rc1=[0; rc1y; 0];
rc2=[rc2x; 0; 0];

syms dc2 real
r0c1=[0; q1; 0];
r0c2=[dc2*cos(q2); q1+ dc2*sin(q2);0];

U1=simplify(-m1*[0, -g0, 0]*r0c1);
U2=simplify(-m2*[0, -g0, 0]*r0c2);
U=simplify(U1+U2);

%the gravity term is:
g_q=jacobian(U, q);
g_q=expand(simplify(transpose(g_q)))

%the gradient of the gravity term is:
dg_q=jacobian(g_q, q);

%we want the upperbound on the norm of this gradient

n_dg_q=norm(dg_q); % = dc2*g0*m2*sin(q2)

%for q2=+-pi/2 the upperbound is tight
%the tight upper bound for this norm then is:
alfa=dc2*g0*m2
