clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 a4 a5 a6 m1 m2 m3 m4 g0 q4 real
syms q1_dot q2_dot q3_dot q4_dot q1_dotdot q2_dotdot q3_dotdot q4_dotdot real
syms rc1x rc1y rc2x rc2y rc3x rc3y % center of mass not on the link axis

M=[a1+2*a2*q2+a3*q2^2+2*a4*q2*sin(q3)+a5*sin(q3)^2, 0, 0;
    0, a3, a4*cos(q3);
    0, a4*cos(q3), a6];

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

c=simplify([c1;c2;c3])


%ii)
M_dot=diff(M,q1)*q1_dot+diff(M,q2)*q2_dot+diff(M,q3)*q3_dot;
S1=q_dot'*C1;
S2=q_dot'*C2;
S3=q_dot'*C3;
S=[S1;S2;S3];

%check if skew symmetric:
k1=M_dot-2*S %it is skew symmetric

%now let's compute another S (S') such that M_dot-2S' is still
%skew-symmetric. If we sum a skew-symmetric matrix to S, the result is still
%skew-symmetric:

S_p=S+[0 -q3_dot q2_dot;q3_dot 0 -q1_dot;-q2_dot q1_dot 0];
k2=M_dot-2*S_p  %it is skew symmetric
%it is a valid choice because [0 -q3_dot q2_dot;q3_dot 0 -q1_dot;-q2_dot
%q1_dot 0]*q_dot=0


S_pp=S+[0 -q3_dot q2_dot;q3_dot 0 -q1_dot; 0 0 0];
k3=M_dot-2*S_pp  %it is NOT skew symmetric

%iii)
q_dotdot=[q1_dotdot; q2_dotdot;q3_dotdot];
u=M*q_dotdot+c; %dynamic model of the robot
a=[a1;a2;a3;a4;a5;a6]; %dynamic coefficients
Y=simplify(jacobian(u,a)) %regressor matrix