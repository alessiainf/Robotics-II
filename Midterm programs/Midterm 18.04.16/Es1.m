clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q4 d4 q3 al3 d3 a3 m1 m2 m3 m4 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot q4_dot real
syms q1_dotdot q2_dotdot q3_dotdot q4_dotdot real
syms Ic1 Ic2 Ic3 Ic4 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real

%prismatic joints
v1=q1_dot;
T1=1/2*m1*norm(v1)^2;

v2=[q1_dot; q2_dot;0];
T2=1/2*m2*norm(v2)^2;

%revolute joints
%direct kinematics
w3=q3_dot;
pc3=[q1+d3*cos(q3);
    q2+d3*sin(q3)];
pc3_d=[q1_dot-d3*sin(q3)*q3_dot;
    q2_dot+d3*cos(q3)*q3_dot]; %vc3
vc3=pc3_d;
T3=1/2*m3*norm(vc3)^2+1/2*w3'*Ic3*w3;

w4=q3_dot+q4_dot;
pc4=[q1+l3*cos(q3)+d4*cos(q3+q4);
    q2+l3*sin(q3)+d4*sin(q3+q4)];
pc4_d=[q1_dot-l3*sin(q3)*q3_dot-d4*sin(q3+q4)*(q3_dot+q4_dot);
    q2_dot+l3*cos(q3)*q3_dot+d4*cos(q3+q4)*(q3_dot+q4_dot)]; %vc4
vc4=pc4_d;
T4=1/2*m4*norm(vc4)^2+1/2*w4'*Ic4*w4;

T=simplify(T1+T2+T3+T4);
T=collect(T, [q1_dot, q2_dot, q3_dot, q4_dot]);


%(m1/2 + m2/2 + m3/2 + m4/2)*q1_dot^2 
% (- (m4*(2*d4*sin(q3 + q4) + 2*l3*sin(q3)))/2 - d3*m3*sin(q3))*q1_dot*q3_dot 
% - d4*m4*sin(q3 + q4)*q1_dot*q4_dot 
% + (m2/2 + m3/2 + m4/2)*q2_dot^2 
% + ((m4*(2*d4*cos(q3 + q4) + 2*l3*cos(q3)))/2 + d3*m3*cos(q3))*q2_dot*q3_dot 
% + d4*m4*cos(q3 + q4)*q2_dot*q4_dot 
% + (Ic3/2 + Ic4/2 + (m4*((d4*cos(q3 + q4) + l3*cos(q3))^2 + (d4*sin(q3 + q4) + l3*sin(q3))^2))/2 + (m3*(d3^2*cos(q3)^2 + d3^2*sin(q3)^2))/2)*q3_dot^2 
% + (Ic4 + (m4*(2*d4*cos(q3 + q4)*(d4*cos(q3 + q4) + l3*cos(q3)) + 2*d4*sin(q3 + q4)*(d4*sin(q3 + q4) + l3*sin(q3))))/2)*q3_dot*q4_dot 
% + (Ic4/2 + (m4*(d4^2*cos(q3 + q4)^2 + d4^2*sin(q3 + q4)^2))/2)*q4_dot^2

m11=(m1 + m2 + m3 + m4);
m12=0;
m13=simplify((- (m4*(d4*sin(q3 + q4) + l3*sin(q3))) - d3*m3*sin(q3)));
m13=collect(m13, [sin(q3+q4), sin(q3)]);
m14=-d4*m4*sin(q3 + q4);

m21=m12;
m22=(m2 + m3 + m4);
m23=simplify(((m4*(d4*cos(q3 + q4) + l3*cos(q3))) + d3*m3*cos(q3)));
m23=collect(m23, [cos(q3+q4), cos(q3)]);
m24=d4*m4*cos(q3 + q4);

m31=m13;
m32=m23;
m33=simplify(Ic3 + Ic4 + (m4*((d4*cos(q3 + q4) + l3*cos(q3))^2 + (d4*sin(q3 + q4) + l3*sin(q3))^2)) + (m3*(d3^2*cos(q3)^2 + d3^2*sin(q3)^2)));
m33=collect(m33, [cos(q3+q4), cos(q3), sin(q3+q4), sin(q3)]);
m34=simplify (Ic4 + (m4*(d4*cos(q3 + q4)*(d4*cos(q3 + q4) + l3*cos(q3)) + d4*sin(q3 + q4)*(d4*sin(q3 + q4) + l3*sin(q3)))));
m34=collect(m34, [cos(q3+q4), cos(q3), sin(q3+q4), sin(q3)]);
m41=m14;
m42=m24;
m43=m34;
m44=simplify(Ic4 + (m4*(d4^2*cos(q3 + q4)^2 + d4^2*sin(q3 + q4)^2)));
m44=collect(m44, [cos(q3 + q4), sin(q3 + q4)]);

M=[m11, m12, m13, m14;
    m21,m22, m23, m24;
    m31, m32, m33, m34;
    m41, m42, m43, m44];

syms a1 a2 a3 a4 a5 a6 real
%i define the dynamic terms
%a1=(m1 + m2 + m3 + m4)
%a2=(d4*m4)
%a3=(d3*m3 + l3*m4)
%a4=m2 + m3 + m4
%a5=m3*d3^2 + m4*d4^2+m4*l3^2 + Ic3 + Ic4
%a6=m4*d4^2+ Ic4

m11=subs(m11, m1 + m2 + m3 + m4, a1);
m13=subs(m13, [-d4*m4, - d3*m3 - l3*m4], [-a2, -a3]);
m14=subs(m14, -d4*m4, -a2);
m22=subs(m22, m2 + m3 + m4, a4);
m23=subs(m23, [d4*m4, d3*m3 + l3*m4],[a2, a3]);
m24=subs(m24, d4*m4, a2);
m31=m13;
m32=m23;
m33=a5+2*a2*l3*cos(q4);
%m33=subs(m33, [ m3*d3^2 + m4*d4^2+ m4*l3^2 + Ic3 + Ic4, 2*m4*cos(q4)*d4*l3], [a5, 2*a2*l3*cos(q4)]);
%m34=subs(m34, [m4*d4^2+Ic4, l3*m4*cos(q4)*d4], [a6, a2*l3*cos(q4)]);
m34=a6+a2*l3*cos(q4);
m41=m14;
m42=m24;
m43=m34;
%m44=subs(m44, m4*d4^2 + Ic4, a6);
m44=a6;

M=[m11, m12, m13, m14;
    m21,m22, m23, m24;
    m31, m32, m33, m34;
    m41, m42, m43, m44]


%centrifugal terms
q=[q1 q2 q3 q4];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));
M3=M(:,3);
C3=(1/2)*(jacobian(M3,q)+jacobian(M3,q)'-diff(M,q3));
M4=M(:,4);
C4=(1/2)*(jacobian(M4,q)+jacobian(M4,q)'-diff(M,q4));

q_dot=[q1_dot; q2_dot; q3_dot; q4_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;
c3=q_dot'*C3*q_dot;
c4=q_dot'*C4*q_dot;

c=simplify([c1;c2;c3; c4])


%gravity terms
syms d1 d2 d3 d4 real
r0c1=[q1+d1; 0; 0];
r0c2=[q1; q2-d2;0];
r0c3=[d3*cos(q3); q2+ d3*sin(q3);0];
r0c4=[dc2*cos(q2); q2+ l3*sin(q3)+d4*sin(q3+q4);0];

U1=simplify(-m1*[0, -g0, 0]*r0c1);
U2=simplify(-m2*[0, -g0, 0]*r0c2);
U3=simplify(-m3*[0, -g0, 0]*r0c3);
U4=simplify(-m4*[0, -g0, 0]*r0c4);
U=simplify(U1+U2+U3+U4);

%the gravity term is:
g_q=jacobian(U, q);
g_q=(simplify(transpose(g_q)));
g_q=[0;
    g0*a4;
    g0*a2*cos(q3+q4)+g0*cos(q3)*a3;
    g0*cos(q3+q4)*a2]

%now i can add the viscous term
syms fv1 fv2 fv3 fv4 a7 a8 a9 a10 real
vis=[fv1*q1_dot; fv2*q2_dot;fv3*q3_dot;fv4*q4_dot];
%defining new dynamic terms:
vis=[a7*q1_dot; a8*q2_dot;a9*q3_dot;a10*q4_dot];

%the dynamic model is
a=[a1; a2; a3; a4; a5; a6; a7; a8; a9; a10];
q_dotdot=[q1_dotdot; q2_dotdot; q3_dotdot; q4_dotdot];
u=M*q_dotdot+c+g_q+vis;

%the regression matrix is:
Y=simplify(jacobian(u,a)) 