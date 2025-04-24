clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 real
syms q1_dot q2_dot q3_dot real
syms rc1x rc2x rc2y rc3x rc3y real% center of mass not on the link axis
syms Ic1xx Ic1yy Ic1zz Ic2xx Ic2yy Ic2zz Ic3xx Ic3yy Ic3zz real

q=[q1, q2, q3];
H_1 = [
    cos(q1), -sin(q1)*cos(al1), sin(q1)*sin(al1), a1*cos(q1);
    sin(q1), cos(q1)*cos(al1), -cos(q1)*sin(al1), a1*sin(q1);
    0, sin(al1), cos(al1), d1;
    0, 0, 0, 1];
H1_ = subs(H_1, [al1, a1, d1, q1], [0, l1, 0, q1]);
R0_1=H1_(1:3, 1:3);
r0_01=H1_(1:3, end);
H_2 = [
    cos(q2), -sin(q2)*cos(al2), sin(q2)*sin(al2), a2*cos(q2);
    sin(q2), cos(q2)*cos(al2), -cos(q2)*sin(al2), a2*sin(q2);
    0, sin(al2), cos(al2), d2;
    0, 0, 0, 1];
H2_ = subs(H_2, [al2, a2, d2, q2], [pi/2, 0, l2, q2]);
R1_2=H2_(1:3, 1:3);
r1_12=H2_(1:3, end);
H_3 = [
    cos(q3), -sin(q3)*cos(al3), sin(q3)*sin(al3), a3*cos(q3);
    sin(q3), cos(q3)*cos(al3), -cos(q3)*sin(al3), a3*sin(q3);
    0, sin(al3), cos(al3), d3;
    0, 0, 0, 1];
H3_ = subs(H_3, [al3, a3, d3, q3], [0, l3, 0, q3]);
R2_3=H3_(1:3, 1:3);
r2_23=H3_(1:3, end);

%kinematic energy of link 1
rc1=[0;rc2y;0];
Ic1=[Ic1xx,0,0;
    0,Ic1yy,0;
    0,0,Ic1zz];
w1=transpose(R0_1)*(q1_dot*[0;0;1]);
v1=transpose(R0_1)*(0+cross(R0_1*w1, r0_01));
vc1=v1+cross(w1,rc1);
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*Ic1*w1);
%T1=collect(T1, q1_dot);

%kinematic energy of link 2
rc2=[0;rc2y;0];
Ic2=[Ic2xx,0,0;
    0,Ic2yy,0;
    0,0,Ic2zz];
w2=transpose(R1_2)*(w1+q2_dot*[0;0;1]);
v2=transpose(R1_2)*(v1++cross(R1_2*w2, r1_12));
vc2=v2+cross(w2,rc2);
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*Ic2*w2);
%T2=collect(T2, [q2_dot, q1_dot]);


%kinematic energy of link 3
rc3=[rc3x;0;0];
Ic3=[Ic3xx,0,0;
    0,Ic3yy,0;
    0,0,Ic3zz];
w3=transpose(R2_3)*(w2+q3_dot*[0;0;1]);
v3=transpose(R2_3)*(v2++cross(R2_3*w3, r2_23));
vc3=v3+cross(w3,rc3);
T3=simplify(1/2*m3*norm(vc3)^2+1/2*w3'*Ic3*w3);
%T3=collect(T3, [q2_dot, q1_dot])

%total kinematic energy
T=T1+T2+T3;
%T=collect(T, [q1_dot, q2_dot, q3_dot]);
T=subs(T, [sin(q2),sin(q3)], [1-cos(q2)^2, 1-cos(q3)^2]);
T=collect(T, [q1_dot, q2_dot, q3_dot]);


syms a1 a2 a3 a4 a5 a6
%the expression of the total kinematic energy is:
% (Ic3xx/2 + Ic2yy/2 + Ic1zz/2 + (l1^2*m1)/2 + (l1^2*m2)/2 + (m1*rc2y^2)/2 + cos(q3)^2*(Ic3yy/2 - Ic3xx/2 + (m3*(l3 + rc3x)^2)/2) + 
% (l1^2*m3*cos(q2)^2)/2 + (l1^2*m3*(cos(q2)^2 - 1)^2*(cos(q3)^2 - 1)^2)/2 + (l1^2*m3*cos(q3)^2*(cos(q2)^2 - 1)^2)/2 + l1*m3*cos(q2)*cos(q3)*(l3 + rc3x))*q1_dot^2 + 
% 
% +(Ic3xx + Ic2yy + cos(q3)^2*(Ic3yy - Ic3xx + m3*(l3 + rc3x)^2) + l1*m3*cos(q2)*cos(q3)*(l3 + rc3x))*q1_dot*q2_dot -
% 
% - l1*m3*(l3 + rc3x)*(cos(q2)^2 - 1)*(cos(q3)^2 - 1)*q1_dot*q3_dot +
% 
% + (Ic3xx/2 + Ic2yy/2 + cos(q3)^2*(Ic3yy/2 - Ic3xx/2 + (m3*(l3 + rc3x)^2)/2))*q2_dot^2 +
% 
% +(Ic3zz/2 + (m3*(l3 + rc3x)^2)/2)*q3_dot^2


m11=(Ic3xx + Ic2yy + Ic1zz + (l1^2*m1) + (l1^2*m2) + (m1*rc2y^2) + cos(q3)^2*(Ic3yy - Ic3xx + (m3*(l3 + rc3x)^2)) + (l1^2*m3*cos(q2)^2) + (l1^2*m3*(cos(q2)^2 - 1)^2*(cos(q3)^2 - 1)^2) + (l1^2*m3*cos(q3)^2*(cos(q2)^2 - 1)^2)/2 + l1*m3*cos(q2)*cos(q3)*(l3 + rc3x));
m11=simplify(m11);

m12=(Ic3xx + Ic2yy + cos(q3)^2*(Ic3yy - Ic3xx + m3*(l3 + rc3x)^2) + l1*m3*cos(q2)*cos(q3)*(l3 + rc3x));
m12=simplify(m12);

m13=- l1*m3*(l3 + rc3x)*(cos(q2)^2 - 1)*(cos(q3)^2 - 1);
m13=simplify(m13);

m22=(Ic3xx + Ic2yy + cos(q3)^2*(Ic3yy - Ic3xx + (m3*(l3 + rc3x)^2)));
m22=simplify(m22);

m23=0;

m33=(Ic3zz + (m3*(l3 + rc3x)^2));
m33=simplify(m33);

%the final expression of the inertia matrix is:
M=[m11, m12, m13;
    m12, m22, m23;
     m13, m23, m33]