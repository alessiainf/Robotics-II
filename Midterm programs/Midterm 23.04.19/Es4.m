clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 real
syms I1 I2 I3 real
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
r0_01=H1_(1:3, end);

H_2 = [
    cos(q2), -sin(q2)*cos(al2), sin(q2)*sin(al2), a2*cos(q2);
    sin(q2), cos(q2)*cos(al2), -cos(q2)*sin(al2), a2*sin(q2);
    0, sin(al2), cos(al2), d2;
    0, 0, 0, 1
];

H2_ = subs(H_2, [al2, a2, d2, q2], [-pi/2, 0, q2, pi/2]);
R1_2=H2_(1:3, 1:3);
r1_12=H2_(1:3, end);

H_3 = [
    cos(q3), -sin(q3)*cos(al3), sin(q3)*sin(al3), a3*cos(q3);
    sin(q3), cos(q3)*cos(al3), -cos(q3)*sin(al3), a3*sin(q3);
    0, sin(al3), cos(al3), d3;
    0, 0, 0, 1
];

H3_ = subs(H_3, [al3, a3, d3, q3], [0, l3, 0, q3]);
R2_3=H3_(1:3, 1:3);
r2_23=H3_(1:3, end);

rc1=[0; rc1y; 0];
rc2=[0; rc2y; 0];
rc3=[rc3x; 0; 0];

w1=transpose(R0_1)*(0+q1_dot*[0;0;1]);
v1=transpose(R0_1)*(0+cross((R0_1*w1), r0_01));
vc1=v1+cross(w1,rc1);
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*I1*w1);

w2=transpose(R1_2)*(w1);
v2=transpose(R1_2)*(v1+q2_dot*[0;0;1]+cross((R1_2*w2), r1_12));
vc2=v2+cross(w2,rc2);
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*I2*w2);

w3=transpose(R2_3)*(w2+q3_dot*[0;0;1]);
v3=transpose(R2_3)*(v2+cross((R2_3*w2), r2_23));
vc3=v3+cross(w3,rc3);
T3=simplify(1/2*m3*norm(vc3)^2+1/2*w3'*I3*w3);

T=T1+T2+T3;
T=collect(T, [q1_dot, q2_dot, q3_dot]);


%inertia matrix
q_dot=[q1_dot; q2_dot; q3_dot];
M=simplify(hessian(T,q_dot))

