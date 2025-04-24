clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 m4 g0 q4 real
syms q1_dot q2_dot q3_dot q4_dot q1_dotdot q2_dotdot q3_dotdot q4_dotdot real
syms rc1x rc1y rc2x rc2y rc3x rc3y % center of mass not on the link axis

vc1=[q1_dot; 0];
T1=simplify(1/2*m1*norm(vc1)^2);

vc2=[q1_dot; q2_dot];
T2=simplify(1/2*m2*norm(vc2)^2);

vc3=[q1_dot+q3_dot; q2_dot];
T3=simplify(1/2*m3*norm(vc3)^2);

vc4=[q1_dot+q3_dot; q2_dot+q4_dot];
T4=simplify(1/2*m4*norm(vc4)^2);

T=simplify(T1+T2+T3+T4);
T=collect(T, [q1_dot, q2_dot, q3_dot, q4_dot]);

%inertia matrix
q_dot=[q1_dot; q2_dot; q3_dot; q4_dot];
M=simplify(hessian(T,q_dot))


%joint velocity
syms vxd vyd real
vd=[vxd;vyd];
p=[q1+q3;q2+q4];
q=[q1;q2;q3;q4];
J=jacobian(p, q);

%The jacobian is full rank but the robot is redundant, so we need to use
%the pseudoinverse instead of the inverse. Also we want to minimize the
%kinematic energy of the robot so we can use M as the weight matrix:
Jm_p=inv(M)*J'*inv((J*inv(M)*J'));
qd_dot=simplify(Jm_p*vd)

% qd_dot =
%   0
%   0
% vxd
% vyd -> only move the last two joints

%the minimum norm solution is given instead by:
qd_dot=simplify(pinv(J)*vd)

% qd_dot =
% vxd/2
% vyd/2
% vxd/2
% vyd/2 -> equally distribute the velocity between joints


