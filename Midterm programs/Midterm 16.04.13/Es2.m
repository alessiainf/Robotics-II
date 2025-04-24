clc; clear all
clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 q4 al2 d2 a2 l3 q3 al3 l4 d3 m4 a3 m1 m2 a1n a2n q1n q2n m3 g0 t1 t2 real
syms q1_dot q2_dot q3_dot real
syms rc1x rc2x rc2y rc3x rc3y real% center of mass not on the link axis
syms Ic1xx Ic1yy Ic1zz Ic2xx Ic2yy Ic2zz Ic3xx Ic3yy Ic3zz real

% direct kinematics of 4R planar robot:
p=[l1*cos(q1)+l2*cos(q1+q2)+l3*cos(q1+q2+q3)+l4*cos(q1+q2+q3+q4);
    l1*sin(q1)+l2*sin(q1+q2)+l3*sin(q1+q2+q3)+l4*sin(q1+q2+q3+q4)];

q=[q1;q2;q3;q4];

%i)
syms d1 d2 d3 d4 real
%only y-components because the other components will be 0
r0c1=d1*sin(q1);
r0c2=l1*sin(q1)+d2*sin(q1+q2);
r0c3=l1*sin(q1)+l2*sin(q1+q2)+d3*sin(q1+q2+q3);
r0c4=l1*sin(q1)+l2*sin(q1+q2)+l3*sin(q1+q2+q3)+d4*sin(q1+q2+q3+q4);

U1=simplify(-m1*-g0*r0c1);
U2=simplify(-m2*-g0*r0c2);
U3=simplify(-m3*-g0*r0c3);
U4=simplify(-m4*-g0*r0c4);
U=simplify(U1+U2+U3+U4);

%the gravity term is:
g_q=jacobian(U, q);
g_q=(simplify(transpose(g_q)));
g_q=collect(g_q, [cos(q1), cos(q1 + q2), cos(q1 + q2 + q3), cos(q1 + q2 + q3 + q4)])

%ii)
%g_q(4)=0 -> cos(q1 + q2 + q3 + q4) =0
%g_q(3)=0 -> cos(q1 + q2 + q3), cos(q1 + q2 + q3 + q4) =0
%g_q(2)=0 -> cos(q1 + q2), cos(q1 + q2 + q3), cos(q1 + q2 + q3 + q4) =0
%g_q(1)=0 -> cos(q1), cos(q1 + q2), cos(q1 + q2 + q3), cos(q1 + q2 + q3 +q4) =0

% -> qe s.t.  cos(q1), cos(q1 + q2), cos(q1 + q2 + q3), cos(q1 + q2 + q3 + q4) =0

%iii)
syms g1 g2 g3 g4 real
%g1=(d1*g0*m1 + g0*l1*m2 + g0*l1*m3 + g0*l1*m4)
%g2=(d2*g0*m2 + g0*l2*m3 + g0*l2*m4)
%g3=(d3*g0*m3 + g0*l3*m4)
%g4=d4*g0*m4

g=[g1;g2;g3;g4];
g_qq=subs(g_q, [(d1*g0*m1 + g0*l1*m2 + g0*l1*m3 + g0*l1*m4), (d2*g0*m2 + g0*l2*m3 + g0*l2*m4), (d3*g0*m3 + g0*l3*m4),  d4*g0*m4 ], [g1, g2, g3, g4]);
Yg=jacobian(g_qq, g)


%iv)
% d4=0;
% d3=-(m4*l3)/m3;
% d2=-(m3+m4)*l2/m2;
% d1=-(m2+m3+m4)*l1/m1;

g_q0=simplify(subs(g_q, [d1, d2, d3, d4], [-((m2+m3+m4)*l1)/m1, -((m3+m4)*l2)/m2, -(m4*l3)/m3, 0]))