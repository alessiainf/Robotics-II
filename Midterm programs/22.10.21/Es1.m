clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real
syms a1 a2 a3 real

%1)
w1=0; 
vc1=[-sin(q1)*q1_dot; 
    cos(q1)*q1_dot]; 
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*I1*w1);

w2=q2_dot;
%i find the velocity of the center of mass as the derivative of the
%position of the center of mass
vc2=[q1_dot-dc2*sin(q2)*q2_dot; 
    dc2*cos(q2)*q2_dot]; 
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*I2*w2);

T=expand(T1+T2);
T=collect(T, [q1_dot, q2_dot]);

q_dot=[q1_dot; q2_dot];
%intertia matrix
M=simplify(hessian(T,q_dot))


% %Now i find the coriolis and centrifugal terms:
q=[q1 q2];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));

q_dot=[q1_dot; q2_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;

c=simplify([c1;c2]) %0
 
%now i find the gravity terms
syms any real
r0c1=[q1 ; any ; any];
r0c2=[q1+dc2*cos(q2);any ; any];

U1=simplify(-m1*[g0, 0, 0]*r0c1);
U2=simplify(-m2*[g0, 0, 0]*r0c2);
U=simplify(U1+U2);

%the gravity term is:
g_q=jacobian(U, q);
g_q=expand(simplify(transpose(g_q)));

%2)
syms a1 a2 a3 real
%a1=m1+m2
%a2=dc2*m2
%a3=m2*dc2^2 + I2

M=subs(M, [m1+m2, m2*dc2^2 + I2,dc2*m2], [a1,a3,a2])
c=subs(c, [m1+m2,dc2*m2, m2*dc2^2 + I2], [a1,a2,a3])
g_q=subs(g_q, [- g0*m1 - g0*m2,dc2*m2, m2*dc2^2 + I2], [-a1*g0,a2,a3])

%Now i can build the regressor matrix as:
q_dotdot=[q1_dotdot; q2_dotdot];
u=M*q_dotdot+c+g_q; %dynamic model of the robot
a=[a1;a2;a3]; %dynamic coefficients
Y=simplify(jacobian(u,a)) %regressor matrix
