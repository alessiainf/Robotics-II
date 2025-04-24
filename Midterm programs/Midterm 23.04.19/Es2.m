clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0
syms w1_1 v1_1 vc1  
syms q1_dot q2_dot q3_dot
syms rc1x rc1y rc2x rc2y rc3x rc3y % center of mass not on the link axis

q=[q1, q2, q3];

H_1 = [
    cos(q1), -sin(q1)*cos(al1), sin(q1)*sin(al1), a1*cos(q1);
    sin(q1), cos(q1)*cos(al1), -cos(q1)*sin(al1), a1*sin(q1);
    0, sin(al1), cos(al1), d1;
    0, 0, 0, 1
];
H1_ = subs(H_1, [al1, a1, d1, q1], [0, l1, 0, q1]);
R0_1=H1_(1:3, 1:3);

H_2 = [
    cos(q2), -sin(q2)*cos(al2), sin(q2)*sin(al2), a2*cos(q2);
    sin(q2), cos(q2)*cos(al2), -cos(q2)*sin(al2), a2*sin(q2);
    0, sin(al2), cos(al2), d2;
    0, 0, 0, 1
];
H2_ = subs(H_2, [al2, a2, d2, q2], [0, l2, 0, q2]);
R1_2=H2_(1:3, 1:3);

H_3 = [
    cos(q3), -sin(q3)*cos(al3), sin(q3)*sin(al3), a3*cos(q3);
    sin(q3), cos(q3)*cos(al3), -cos(q3)*sin(al3), a3*sin(q3);
    0, sin(al3), cos(al3), d3;
    0, 0, 0, 1
];
H3_ = subs(H_3, [al3, a3, d3, q3], [0, l3, 0, q3]);
R2_3=H3_(1:3, 1:3);

rc1=[rc1x; rc1y; 0];
rc2=[rc2x; rc2y; 0];
rc3=[rc3x; rc3y; 0];


%now i find the gravity terms
syms any real
r0c1=H1_*[rc1;1]; 
r0c1=r0c1(1:3);
r0c2=H1_*H2_*[rc2;1];
r0c2=r0c2(1:3);
r0c3=H1_*H2_*H3_*[rc3;1];
r0c3=r0c3(1:3);

U1=simplify(-m1*[g0, 0, 0]*r0c1);
U2=simplify(-m2*[g0, 0, 0]*r0c2);
U3=simplify(-m3*[g0, 0, 0]*r0c3);
U=simplify(U1+U2+U3);

%the gravity term is:
g_q=jacobian(U, q);
g_q=(simplify(transpose(g_q)));

%in order to have the same g(q) requested, we can first substitute some
%values to have 0 in the second and third rows:

g_q=subs(g_q, [rc3x, rc3y], [-l3, 0]); %third row
g_q=simplify(subs(g_q, [rc2y, rc2x], [0, (-(m2+m3)*l2)/m2]));  %second row
g_q=simplify(subs(g_q, rc1x, -(m1+m2+m3)*l1/m1)) %first row

















