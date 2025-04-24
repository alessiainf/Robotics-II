clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real
syms a1 a2 a3 real

%the robot is moving in an horizontal plane -> g_q=0

M = [a2*q2^2 + a1, -a2*l1;
     -a2*l1,     a2];
 
 
c = [2*a2*q2*q1_dot*q2_dot;
      -a2*q2*q1_dot^2];

a=[a1;a2;a3]; %dynamic coefficients
% a1= m1*dc1^2 + m2*l1^2+ I1 + I2;
% a2=m2;
% a3=dc1*m1+l1*m2;

%dynamic model of the robot is:
q_dotdot=[q1_dotdot; q2_dotdot];
u=M*q_dotdot+c;
%regressor matrix
%Y=simplify(jacobian(u,a)) %regressor matrix

%the trajectory is:
syms a b k t T
qd=[a+b*(1-cos((pi*t)/T)); k];
qd_dot=diff(qd, t);
qd_dotdot=diff(qd_dot, t);

%the torque is:
tau=subs(u, [q1, q2, q1_dot, q2_dot, q1_dotdot, q2_dotdot], [qd(1), qd(2), qd_dot(1), qd_dot(2), qd_dotdot(1), qd_dotdot(2)]);
%tau=simplify(subs(tau, [a1,a2,a3], [m1*dc1^2 + m2*l1^2+ I1 + I2, m2, dc1*m1+l1*m2]))

%the maximum absolute value of the torque is reached at the initial and
%final instants -> t=0, t=T
syms tau_max real
tau1=subs(tau(1), t, 0) %=-tau_1=subs(tau(1), t, T) -> tau1(1) need to be lower than the bound t1_max
% tau1 =
% (b*pi^2*(a2*k^2 + a1))/T^2
%we found the minimum time T by isolating T in the first term of tau1
%because of the bound
T1=sqrt((b*pi^2*(a2*k^2 + a1))/tau_max)
%the first joint saturates at time 0 and T, so the value of the secondo
%torque is:
tau2=subs(tau(2), t, 0) %=-tau_1=subs(tau(2), t, T) 