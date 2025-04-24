clc; clear all
clc; clear all;
syms l q1 al1 d1 a1 l2 q2 q4 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 a1n a2n q1n q2n m3 g0 t1 t2 real
syms q1_dot q2_dot q3_dot real
syms rc1x rc2x rc2y rc3x rc3y real% center of mass not on the link axis
syms I1 I2 I3 real

p=[q1+l2*cos(q2)+l3*cos(q2+q3);
    l2*sin(q2)+l3*sin(q2+q3)];

w1=0; %prismatic joint
vc1=q1_dot; 
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*I1*w1);

w2=q2_dot;
%i find the velocity of the center of mass as the derivative of the
%position of the center of mass
vc2=[q1_dot-d2*sin(q2)*q2_dot;
    d2*cos(q2)*q2_dot];
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*I2*w2);

w3=q2_dot+q3_dot;
%w3=w2; %prismatic joint -> depends on the previous joint
%i find the velocity of the center of mass as the derivative of the
%position of the center of mass
vc3=[q1_dot-l2*sin(q2)*q2_dot-d3*sin(q2+q3)*(q2_dot+q3_dot);    
     l2*cos(q2)*q2_dot+d3*cos(q2+q3)*(q2_dot+q3_dot)];
T3=simplify(1/2*m3*norm(vc3)^2+1/2*w3'*I3*w3);

T=expand(T1+T2+T3);
T=collect(T, [q1_dot, q2_dot, q3_dot]);

%total kinematic energy
% ((m1 + m2 + m3)*q1_dot^2)/2 
% + (- d2*m2*sin(q2) - l2*m3*sin(q2) - d3*m3*cos(q2)*sin(q3) - d3*m3*cos(q3)*sin(q2))*q1_dot*q2_dot 
% + (- d3*m3*cos(q2)*sin(q3) - d3*m3*cos(q3)*sin(q2))*q1_dot*q3_dot 
%((m2*d2^2) + (m3*d3^2*cos(q2)^2*cos(q3)^2) + (m3*d3^2*cos(q2)^2*sin(q3)^2) + (m3*d3^2*cos(q3)^2*sin(q2)^2) + (m3*d3^2*sin(q2)^2*sin(q3)^2) + m3*d3*l2*cos(q2)^2*cos(q3) + m3*d3*l2*cos(q3)*sin(q2)^2 + (m3*l2^2*cos(q2)^2) + (m3*l2^2*sin(q2)^2) + I2 + I3)*q2_dot^2  
% + (m3*d3^2*cos(q2)^2*cos(q3)^2 + m3*d3^2*cos(q2)^2*sin(q3)^2 + m3*d3^2*cos(q3)^2*sin(q2)^2 + m3*d3^2*sin(q2)^2*sin(q3)^2 + l2*m3*d3*cos(q2)^2*cos(q3) + l2*m3*d3*cos(q3)*sin(q2)^2 + I3)*q2_dot*q3_dot 
% + ((m3*d3^2*cos(q2)^2*cos(q3)^2 + m3*d3^2*cos(q2)^2*sin(q3)^2 + m3*d3^2*cos(q3)^2*sin(q2)^2 + m3*d3^2*sin(q2)^2*sin(q3)^2 + I3)*q3_dot^2)/2


m11=simplify((m1 + m2 + m3));
m12=simplify(- d2*m2*sin(q2) - l2*m3*sin(q2) - d3*m3*cos(q2)*sin(q3) - d3*m3*cos(q3)*sin(q2));
m13=simplify(- d3*m3*cos(q2)*sin(q3) - d3*m3*cos(q3)*sin(q2));
m21=m12;
m22=2*(simplify( ((m2*d2^2)/2 + (m3*d3^2*cos(q2)^2*cos(q3)^2)/2 + (m3*d3^2*cos(q2)^2*sin(q3)^2)/2 + (m3*d3^2*cos(q3)^2*sin(q2)^2)/2 + (m3*d3^2*sin(q2)^2*sin(q3)^2)/2 + m3*d3*l2*cos(q2)^2*cos(q3) + m3*d3*l2*cos(q3)*sin(q2)^2 + (m3*l2^2*cos(q2)^2)/2 + (m3*l2^2*sin(q2)^2)/2 + I2/2 + I3/2)));
m23=simplify( m3*d3^2*cos(q2)^2*cos(q3)^2 + m3*d3^2*cos(q2)^2*sin(q3)^2 + m3*d3^2*cos(q3)^2*sin(q2)^2 + m3*d3^2*sin(q2)^2*sin(q3)^2 + l2*m3*d3*cos(q2)^2*cos(q3) + l2*m3*d3*cos(q3)*sin(q2)^2 + I3);
m31=m13;
m32=m23;
m33=simplify((m3*d3^2*cos(q2)^2*cos(q3)^2 + m3*d3^2*cos(q2)^2*sin(q3)^2 + m3*d3^2*cos(q3)^2*sin(q2)^2 + m3*d3^2*sin(q2)^2*sin(q3)^2 + I3));
M=[m11, m12, m13;
    m21, m22, m23;
    m31, m32, m33];

%dynamic terms
syms a1 a2 a3 a4 a5 a6 real
% a1= m1 + m2 + m3;
% a2=d2*m2+l2*m3;
% a3=m2*d2^2 + I2 +m3*l2^2+m3*d3^2  + I3
% a4= m3*d3^2  + I3
% a5=m3*d3*l2
% a6=d3*m3

m11=a1;
m12=-sin(q2)*a2-a6*sin(q2+q3);
m13=-a6*sin(q2+q3);
m21=m12;
m22= a3+2*a5*cos(q3);
m23=a4+a5*cos(q3);
m31=m13;
m32=m23;
m33=a4;

%intertia matrix
M=[m11, m12, m13;
    m21, m22, m23;
    m31, m32, m33]

%coriolis/centrifugal terms
q=[q1; q2; q3];
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

c=simplify([c1;c2;c3])