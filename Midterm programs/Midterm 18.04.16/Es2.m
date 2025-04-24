clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q4 d4 q3 al3 d3 a3 m1 m2 m3 m4 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot q4_dot real
syms q1_dotdot q2_dotdot q3_dotdot q4_dotdot real
syms Ic1 Ic2 Ic3 Ic4 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real


M=[a1+2*a2*cos(q2), a3+a2*cos(q2);
    a3+a2*cos(q2), a3];

%1)
q=[q1 q2];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));

q_dot=[q1_dot; q2_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;

c=simplify([c1;c2]);

%Now i find the first matrix S
M_dot=diff(M,q1)*q1_dot+diff(M,q2)*q2_dot
S_1=q_dot'*C1;
S_2=q_dot'*C2;
S1=[S_1;S_2];
%check if skew symmetric:
k1=simplify(M_dot-2*S1) %it is skew symmetric

%now let's compute another S (S2) such that M_dot-2S2 is not skew-symmetric. 
% In general, if we sum a skew-symmetric matrix to S1, the result is still
% skew-symmetric. So if we sum a not skew-Symmetric matrix to S1 but such
% that S2*q_dot=0, we can find the requested matrix.

S2=[-2*a2*q2_dot*sin(q2), -a2*q2_dot*sin(q2);
    a2*q2_dot*sin(q2), 0];

k2=simplify(M_dot-2*S2) %it is not skew symmetric


%2)
syms delta_q2 q0_2 q0_1 t T real
%cubic trajectory rest-to-rest:
qd=[q0_1;q0_2]+[0; delta_q2]*(-2*(t/T)^3+3*(t/T)^2); %only the second joint should move -> delta_q=[0; delta_q2]
qd_d=([0; delta_q2]/T)*(-6*(t/T)^2+6*(t/T)); 
qd_dd=([0; delta_q2]/T^2)*(-12*(t/T)+6);

%move in horizontal plane -> g=0
u=M*qd_dd+c

%start at rest -> q_dot=0, t=0
%the initial torque is:
u_0=eval(subs(u, [a1, a2, a3, q2, delta_q2, T, q1_dot, q2_dot, t], [17, 5, 3, -pi/2, pi/2, 2, 0, 0, 0]))

% u_0 =
% 
%     7.0686
%     7.0686

%The second torque move the second joint from -pi/2 to 0 (delta_q2=pi/2) -> positive sign, counterclockwise. 

%on the other hand since the first joint should not move, the first torque is different 
% from zero because it should balance the movement of the second joint due
% to inertial coupling, otherwise it would start move on the opposite direction of the second joint 
% -> there is a positive torque to balance this effect

%also we can see that tau1 depends on (a3 + a2*cos(q2)) -> it depends on the second joints. 
% Since a3 and a2 are positive this influence also the sign of tau1.
%tau2 instead depends only on a3 that it's also positive
