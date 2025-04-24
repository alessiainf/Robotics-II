clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real

w1=q1_dot; %prismatic joint
vc1=[-dc1*sin(q1); dc1*cos(q1)]; 
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*I1*w1);

w2=q1_dot;
%i find the velocity of the center of mass as the derivative of the
%position of the center of mass
vc2=[-l1*sin(q1)*q1_dot+q2_dot*sin(q1)+q2*cos(q1)*q1_dot;
    l1*cos(q1)*q1_dot-q2_dot*cos(q1)+q2*sin(q1)*q1_dot];
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*I2*w2);

T=expand(T1+T2);
T=collect(T, [q1_dot, q2_dot]);

%(m1*dc1^2)/2 + ((m2*l1^2*cos(q1)^2)/2 + (m2*l1^2*sin(q1)^2)/2 + (m2*q2^2*cos(q1)^2)/2 + (m2*q2^2*sin(q1)^2)/2 + I1/2 + I2/2)*q1_dot^2 
% + (- l1*m2*cos(q1)^2 - l1*m2*sin(q1)^2)*q1_dot*q2_dot 
% + ((m2*cos(q1)^2)/2 + (m2*sin(q1)^2)/2)*q2_dot^2




m11=2*simplify((m1*dc1^2)/2 + ((m2*l1^2*cos(q1)^2)/2 + (m2*l1^2*sin(q1)^2)/2 + (m2*q2^2*cos(q1)^2)/2 + (m2*q2^2*sin(q1)^2)/2 + I1/2 + I2/2));

m12=simplify(- l1*m2*cos(q1)^2 - l1*m2*sin(q1)^2);

m21=m12;

m22=2*simplify((m2*cos(q1)^2)/2 + (m2*sin(q1)^2)/2);


%The inertia matrix is:
M=[m11, m12;
    m21, m22]


% %Now i find the coriolis and centrifugal terms:
q=[q1 q2];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));

q_dot=[q1_dot; q2_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;

c=simplify([c1;c2])
 
%now i find the gravity terms
syms any real
r0c1=[any ; dc1*sin(q1); any];
r0c2=[any ; l1*sin(q1)-q2*cos(q1); any];

U1=simplify(-m1*[0, -g0, 0]*r0c1);
U2=simplify(-m2*[0, -g0, 0]*r0c2);
U=simplify(U1+U2);

%the gravity term is:
g_q=jacobian(U, q);
g_q=expand(simplify(transpose(g_q)))


%Now i define the dynamic terms by assuming to kno wl1 and g0:
syms a1 a2 a3 real
% a1= m1*dc1^2 + m2*l1^2+ I1 + I2;
% a2=m2;
% a3=dc1*m1+l1*m2;

m11=subs(m11, [m1*dc1^2 + m2*l1^2 + I1 + I2, m2], [a1, a2]);
m12=subs(m12, [-l1*m2], [-l1*a2]);
m21=m12;
m22=subs(m22, [ m2], [ a2]);
%The new inertia matrix is:
M=[m11, m12;
    m21, m22]

%the new centrifugal term is:
c=subs(c, [m2], [a2])

%the new gravity term is:
g_q=subs(g_q, [dc1*g0*m1*cos(q1) + g0*l1*m2*cos(q1), m2], [a3*g0*cos(q1), a2])

%Now i can build the regressor matrix as:
q_dotdot=[q1_dotdot; q2_dotdot];
u=M*q_dotdot+c+g_q; %dynamic model of the robot
a=[a1;a2;a3]; %dynamic coefficients
Y=simplify(jacobian(u,a)) %regressor matrix


%r is minimal. if dc1=0 however a3 is no longer needed and we only need
%only 2 dynamic terms 


