clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real


M=[a1, a2*sin(q2);
    a2*sin(q2), a3];

%a)

%I find the coriolis and centrifugal terms
q=[q1; q2];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));
q_dot=[q1_dot; q2_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;
c=simplify([c1;c2])

%i check if skew-symmetric
M_dot=diff(M,q1)*q1_dot+diff(M,q2)*q2_dot;
S1=q_dot'*C1
S2=q_dot'*C2
S=[S1;S2];

%check if skew symmetric:
k1=M_dot-2*S %it is skew symmetric

%Now let's find another factorization S_p such that it is not
%skew-simmetric
S_p=[q2_dot -q1_dot;
    0,0];  %it is valid because S_p*q_dot=0
k2=M_dot-2*S_p; % -> it is not skew-symmetric


%b)
% it's possible to always find another factorization matrix S' that
% satisfies S'(q, q_dot)=c(q, q_dot) because S' can be construct as S'=S+S''
% where S'' is a skew-symmetric matrix such that S''*q_dot=0 (i.e. it do not alter the coriolis term). 
% This ensure that M_dot-2*S' remains skew-symmetric.


