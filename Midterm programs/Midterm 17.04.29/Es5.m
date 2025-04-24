clc; clear all
clc; clear all;
syms L1 q1 al1 d1 a1 L2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 a1n a2n q1n q2n m3 g0 t1 t2 real
syms q1_dot q2_dot q3_dot real
syms rc1x rc2x rc2y rc3x rc3y real% center of mass not on the link axis
syms Ic1xx Ic1yy Ic1zz Ic2xx Ic2yy Ic2zz Ic3xx Ic3yy Ic3zz real

%direct kinematics of 2R planar robot:
% r_nom=[L1*cos(q1)+L2*cos(q1+q2);
%     L1*sin(q1)+L2*sin(q1+q2);
%     0];


H_1 = [
    cos(t1), -sin(t1)*cos(al1), sin(t1)*sin(al1), a1*cos(t1);
    sin(t1), cos(t1)*cos(al1), -cos(t1)*sin(al1), a1*sin(t1);
    0, sin(al1), cos(al1), d1;
    0, 0, 0, 1];
H1_ = subs(H_1, [al1, a1, d1, t1], [0, a1, 0, q1]);
H_2 = [
    cos(t2), -sin(t2)*cos(al2), sin(t2)*sin(al2), a2*cos(t2);
    sin(t2), cos(t2)*cos(al2), -cos(t2)*sin(al2), a2*sin(t2);
    0, sin(al2), cos(al2), d2;
    0, 0, 0, 1];
H2_ = subs(H_2, [al2, a2, d2, t2], [0, a2, 0, q2]);

%computed direct kinematics:
r_hom=H1_*H2_*[0;0;0;1];
r=simplify(r_hom(1:2, end));
%->
% r =[L2*cos(q1 + q2) + L1*cos(q1);
% L2*sin(q1 + q2) + L1*sin(q1)];

%nominal direct kinematics:
r_nom=subs(r, [a1, a2, q1, q2], [a1n, a2n, q1n, q2n]);
%->
% r_nom = [a1*cos(t1) + a2*cos(t1)*cos(t2) - a2*sin(t1)*sin(t2);
% a1*sin(t1) + a2*cos(t1)*sin(t2) + a2*cos(t2)*sin(t1)];

%parameters
an=[a1n;a2n];
delta_a=[a1n-a1; a2n-a2];
qn=[q1n, q2n];
delta_q=[q1n-q1;q2n-q2];

%driect kinematics error
%delta_r=r-r_nom
delta_r=jacobian(r_nom, an)*delta_a+jacobian(r_nom, qn)*delta_q;
%-> delta_r=Y*delta_fi

delta_fi=[delta_a; delta_q];
Y=[jacobian(r_nom, an), jacobian(r_nom, qn)] %regressor matrix

%delta_rr=Y*delta_fi %check


