clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real
syms a1 a2 a3 real

%position of center of mass 1
pc1=[(l1+rc1x)*cos(q1); rc1y*sin(q1)];

%position of center of mass 2
pc2=[l1*cos(q1)+(l2+rc2x)*cos(q1+q2)-(rc2y)*sin(q1+q2); 
    l1*sin(q1)+(l2+rc2y)*sin(q1+q2)+(rc2y)*cos(q1+q2)];

w1=q1_dot; 
vc1=[-(l1+rc1x)*sin(q1)*q1_dot; rc1y*cos(q1)*q1_dot];
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*I1*w1);

w2=q1_dot+q2_dot;
%i find the velocity of the center of mass as the derivative of the
%position of the center of mass
vc2=[-l1*sin(q1)*q1_dot+-(l2+rc2x)*sin(q1+q2)*(q1_dot+q2_dot)-(rc2y)*cos(q1+q2)*(q1_dot+q2_dot); 
    l1*cos(q1)*q1_dot+(l2+rc2x)*cos(q1+q2)*(q1_dot+q2_dot)-(rc2y)*sin(q1+q2)*(q1_dot+q2_dot)];

T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*I2*w2);

T=expand(T1+T2);
T=collect(T, [q1_dot, q2_dot])

%inertia matrix
q_dot=[q1_dot; q2_dot];
M=simplify(hessian(T,q_dot));

% I1 + I2 + l1^2*m1 + l1^2*m2 + l2^2*m2 + m1*rc1x^2 + m2*rc2x^2 + m2*rc2y^2 + 2*l1*m1*rc1x + 2*l2*m2*rc2x - l1^2*m1*cos(q1)^2 - m1*rc1x^2*cos(q1)^2 + m1*rc1y^2*cos(q1)^2 + 2*l1*l2*m2*cos(q2) + 2*l1*m2*rc2x*cos(q2) - 2*l1*m2*rc2y*sin(q2) - 2*l1*m1*rc1x*cos(q1)^2, 
% m2*l2^2 + 2*m2*l2*rc2x + l1*m2*cos(q2)*l2 + m2*rc2x^2 + l1*m2*cos(q2)*rc2x + m2*rc2y^2 - l1*m2*sin(q2)*rc2y + I2]
% m2*l2^2 + 2*m2*l2*rc2x + l1*m2*cos(q2)*l2 + m2*rc2x^2 + l1*m2*cos(q2)*rc2x + m2*rc2y^2 - l1*m2*sin(q2)*rc2y + I2,                                                              
% m2*l2^2 + 2*m2*l2*rc2x + m2*rc2x^2 + m2*rc2y^2 + I2

%dynamic terms
% a1=I1 + I2 + l1^2*m1 + l1^2*m2 + l2^2*m2 + m1*rc1x^2 + m2*rc2x^2 + m2*rc2y^2 + 2*l1*m1*rc1x + 2*l2*m2*rc2x
% a2=l1^2*m1- m1*rc1x^2 +m1*rc1y^2- 2*l1*m1*rc1x
% a3=l1*l2*m2+l1*m2*rc2x
% a4=-l1*m2*rc2y
% a5=I2+ m2*rc2y^2 + m2*rc2x^2 +m2*l2^2+  2*m2*l2*rc2x

syms a1 a2 a3 a4 a5 real
m11= a1 - a2*cos(q1)^2  + 2*a3*cos(q2) + 2*a4*sin(q2) ;
m12= a5 + a3*cos(q2) +a4*sin(q2);  
m21= m12;                                                              
m22= a5 ;

M=[m11,m12;
   m21, m22]


% %Now i find the coriolis and centrifugal terms:
q=[q1 q2];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));

c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;

c=simplify([c1;c2]) 
 
%now i find the gravity terms
syms any real
% r0c1=[q1 ; any ; any];
% r0c2=[q1+dc2*cos(q2);any ; any];
r0c1=pc1;
r0c2=pc2;

U1=simplify(-m1*[0, -g0]*r0c1);
U2=simplify(-m2*[0, -g0]*r0c2);
U=simplify(U1+U2);

%the gravity term is:
g_q=jacobian(U, q);
g_q=(simplify(transpose(g_q)));
g_q=collect(g_q, [cos(q1+q2), sin(q1+q2), cos(q1)]);
%dynamic terms
syms a6 a7 a8 a9 real
% a6=g0*m2*(l2 + rc2y);
% a7=(-g0*m2*rc2y);
% a8=(g0*l1*m2 + g0*m1*rc1y);
g_q=[a6*cos(q1 + q2) + a7*sin(q1 + q2) + a8*cos(q1);
                                  a6*cos(q1 + q2) + a7*sin(q1 + q2)]


% %Now i can build the regressor matrix as:
 q_dotdot=[q1_dotdot; q2_dotdot];
 u=M*q_dotdot+c+g_q; %dynamic model of the robot
 a=[a1;a2;a3;a4;a5;a6;a7;a8]; %dynamic coefficients
 Y=simplify(jacobian(u,a)) %regressor matrix