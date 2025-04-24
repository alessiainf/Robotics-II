clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real

w1=0; %prismatic joint
vc1=q1_dot; 
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*I1*w1);

w2=q2_dot;
%i find the velocity of the center of mass as the derivative of the
%position of the center of mass
vc2=[-dc2*sin(q2)*q2_dot;
    q1_dot+dc2*cos(q2)*q2_dot];
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*I2*w2);

w3=w2; %prismatic joint -> depends on the previous joint
%i find the velocity of the center of mass as the derivative of the
%position of the center of mass
vc3=[-(l2*sin(q2)+q3*cos(q2))*q2_dot-sin(q2)*q3_dot;
    q1_dot+(l2*cos(q2)-q3*sin(q2))*q2_dot+cos(q2)*q3_dot];
T3=simplify(1/2*m3*norm(vc3)^2+1/2*w3'*I3*w3);

T=expand(T1+T2+T3);
T=collect(T, [q1_dot, q2_dot, q3_dot]);

% (m1 + m2 + m3)*q1_dot^2 
% + ((m3*(l2*cos(q2) - q3*sin(q2))) + dc2*m2*cos(q2))*q1_dot*q2_dot 
% + m3*cos(q2)*q1_dot*q3_dot 
% + (I2 + I3 + (dc2^2*m2) + (m3*(q3*cos(q2) + l2*sin(q2))^2) + (m3*(l2*cos(q2) - q3*sin(q2))^2))*q2_dot^2 
% + (m3*cos(q2)*(l2*cos(q2) - q3*sin(q2)) + m3*sin(q2)*(q3*cos(q2) + l2*sin(q2)))*q2_dot*q3_dot 
% + ((m3*cos(q2)^2) + (m3*sin(q2)^2))*q3_dot^2



m11=simplify((m1 + m2 + m3));

m12=simplify((m3*(l2*cos(q2) - q3*sin(q2))) + dc2*m2*cos(q2));

m13=simplify(m3*cos(q2));

m21=m12;

m22=simplify(I2 + I3 + (dc2^2*m2) + (m3*(q3*cos(q2) + l2*sin(q2))^2) + (m3*(l2*cos(q2) - q3*sin(q2))^2));

m23=simplify(m3*cos(q2)*(l2*cos(q2) - q3*sin(q2)) + m3*sin(q2)*(q3*cos(q2) + l2*sin(q2)));

m31=m13;

m32=m23;

m33=simplify((m3*cos(q2)^2) + (m3*sin(q2)^2));

%The inertia matrix is:
M=[m11, m12, m13;
    m21, m22, m23;
    m31, m32, m33];

syms a1 a2 a3 a4 real
%Now i define the dynamic terms:
% a1= m1 + m2 + m3;
% a2=m3;
% a3=dc2*m2;
% a4=m2*dc2^2+I2+I3;


%So the new inertia matrix is:
m11=subs(m11, [m1 + m2 + m3], [a1]);
m12=subs(m12, [m3, dc2*m2], [a2, a3]);
m13=subs(m13, [m3,  dc2*m2, I2+I3+m2*dc2^2], [ a2, a3, a4]);
m21=subs(m21, [ m3, dc2*m2, I2+I3+m2*dc2^2], [ a2, a3, a4]);
m22=subs(m22, [ m3,  I2+I3+m2*dc2^2], [ a2,  a4]);
m23=subs(m23, [ m3, dc2*m2, I2+I3+m2*dc2^2], [ a2, a3, a4]);
m31=subs(m31, [ m3, dc2*m2, I2+I3+m2*dc2^2], [ a2, a3, a4]);
m32=subs(m32, [ m3, dc2*m2, I2+I3+m2*dc2^2], [ a2, a3, a4]);
m33=subs(m33, [m3, dc2*m2, I2+I3+m2*dc2^2], [ a2, a3, a4]);

%The new inertia matrix is:
M=[m11, m12, m13;
    m21, m22, m23;
    m31, m32, m33];


%Now i find the coriolis and centrifugal terms:


%i)
q=[q1 q2 q3];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));
M3=M(:,3);
C3=(1/2)*(jacobian(M3,q)+jacobian(M3,q)'-diff(M,q3));

q_dot=[q1_dot; q2_dot; q3_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;
c3=q_dot'*C3*q_dot;

c=simplify([c1;c2;c3]);


%Now i can build the regressor matrix as:
q_dotdot=[q1_dotdot; q2_dotdot; q3_dotdot];
u=M*q_dotdot+c; %dynamic model of the robot
a=[a1;a2;a3;a4]; %dynamic coefficients
Y=simplify(jacobian(u,a)) %regressor matrix


