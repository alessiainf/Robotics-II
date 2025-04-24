clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 A C D E F g0 real
syms Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 Ixx3 Iyy3 Izz3 real
syms q1_dot q2_dot q3_dot real
syms rc1x rc1y rc2x rc2y rc3x rc3y real % center of mass not on the link axis

q=[q1, q2, q3];

H_1 = [
    cos(q1), -sin(q1)*cos(al1), sin(q1)*sin(al1), a1*cos(q1);
    sin(q1), cos(q1)*cos(al1), -cos(q1)*sin(al1), a1*sin(q1);
    0, sin(al1), cos(al1), d1;
    0, 0, 0, 1
];

H1_ = subs(H_1, [al1, a1, d1, q1], [pi/2, 0, l1, q1]);
R0_1=H1_(1:3, 1:3);
r0_01=H1_(1:3, end)

H_2 = [
    cos(q2), -sin(q2)*cos(al2), sin(q2)*sin(al2), a2*cos(q2);
    sin(q2), cos(q2)*cos(al2), -cos(q2)*sin(al2), a2*sin(q2);
    0, sin(al2), cos(al2), d2;
    0, 0, 0, 1
];

H2_ = subs(H_2, [al2, a2, d2, q2], [0, l2, 0, q2]);
R1_2=H2_(1:3, 1:3);
r1_12=H2_(1:3, end)

H_3 = [
    cos(q3), -sin(q3)*cos(al3), sin(q3)*sin(al3), a3*cos(q3);
    sin(q3), cos(q3)*cos(al3), -cos(q3)*sin(al3), a3*sin(q3);
    0, sin(al3), cos(al3), d3;
    0, 0, 0, 1
];

H3_ = subs(H_3, [al3, a3, d3, q3], [0, l3, 0, q3]);
R2_3=H3_(1:3, 1:3);
r2_23=H3_(1:3, end)

rc1=[A; -F; 0];
rc2=[-C; 0; 0];
rc3=[-D; 0; E];

Ic1=[Ixx1, 0 0;
    0, Iyy1,0;
    0,0,Izz1];
Ic2=[Ixx2, 0 0;
    0, Iyy2,0;
    0,0,Izz2];
Ic3=[Ixx3, 0 0;
    0, Iyy3,0;
    0,0,Izz3];

w1=transpose(R0_1)*(0+q1_dot*[0;0;1]);
v1=transpose(R0_1)*(0+cross((R0_1*w1), r0_01));
vc1=v1+cross(w1,rc1);
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*Ic1*w1);

w2=transpose(R1_2)*(w1+q2_dot*[0;0;1]);
v2=transpose(R1_2)*(v1+cross((R1_2*w2), r1_12));
vc2=v2+cross(w2,rc2);
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*Ic2*w2);

w3=transpose(R2_3)*(w2+q3_dot*[0;0;1]);
v3=transpose(R2_3)*(v2+cross((R2_3*w2), r2_23));
vc3=v3+cross(w3,rc3);
T3=simplify(1/2*m3*norm(vc3)^2+1/2*w3'*Ic3*w3);

T=T1+T2+T3;
T=collect(T, [q1_dot, q2_dot, q3_dot]);

%(Ixx2/2 + Ixx3/4 + Iyy1/2 + Iyy3/4 - (Ixx3*cos(2*q2 + 2*q3))/4 + (Iyy3*cos(2*q2 + 2*q3))/4 + (A^2*m1)/2 + (m3*(l2*cos(q2) + l3*cos(q2) - cos(q2 + q3)*D)^2)/2 - (Ixx2*cos(q2)^2)/2 + (Iyy2*cos(q2)^2)/2 + (m2*cos(q2)^2*(C - l2)^2)/2 + (m3*cos(q2 + q3)^2*E^2)/2 + (m3*sin(q2 + q3)^2*E^2)/2)*q1_dot^2 
% + (l2*m3*cos(q2 + q3)*sin(q3)*E - m3*sin(q2 + q3)*E*(l3 - D + l2*cos(q3)))*q1_dot*q2_dot 
% + m3*sin(q2 + q3)*D*E*q1_dot*q3_dot 
% + (Izz2/2 + Izz3/2 + (m2*(C - l2)^2)/2 + (m3*(l3 - D + l2*cos(q3))^2)/2 + (l2^2*m3*sin(q3)^2)/2)*q2_dot^2 
% + (Izz3 - m3*D*(l3 - D + l2*cos(q3)))*q2_dot*q3_dot 
% + (Izz3/2 + (m3*D^2)/2)*q3_dot^2

m11=simplify(Ixx2 + Ixx3/2 + Iyy1 + Iyy3/2 - (Ixx3*cos(2*q2 + 2*q3))/2 + (Iyy3*cos(2*q2 + 2*q3))/2 + (A^2*m1) + (m3*(l2*cos(q2) + l3*cos(q2) - cos(q2 + q3)*D)^2) - (Ixx2*cos(q2)^2) + (Iyy2*cos(q2)^2) + (m2*cos(q2)^2*(C - l2)^2) + (m3*cos(q2 + q3)^2*E^2) + (m3*(1-cos(q2+q3)^2)*E^2));
m11=collect(m11, [cos(q2), cos(q3), cos(q2 + q3)]);
m12=simplify((l2*m3*cos(q2 + q3)*sin(q3)*E - m3*sin(q2 + q3)*E*(l3 - D + l2*cos(q3))));
m12=collect(m12, [cos(q3), sin(q3), cos(q2 + q3), sin(q2 + q3)]);
m13=m3*sin(q2 + q3)*D*E;
m21=m12;
m22=simplify((Izz2 + Izz3 + (m2*(C - l2)^2) + (m3*(l3 - D + l2*cos(q3))^2) + (l2^2*m3*(1-cos(q3)^2))));
m22=collect(m22, [cos(q3)]);
m23=simplify(Izz3 - m3*D*(l3 - D + l2*cos(q3)));
m23=collect(m23, [cos(q3)]);
m31=m13;
m32=m23;
m33=simplify((Izz3 + (m3*D^2)));

%the final intertia matrix:
M=[m11, m12, m13;
   m21, m22, m23;
   m31, m32, m33]