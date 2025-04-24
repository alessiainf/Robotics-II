clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real

%direct kinematics
p=[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
    sin(q1)+ sin(q1+q2)+sin(q1+q2+q3)];
q=[q1 q2 q3];
J=jacobian(p, q);

q_dot=[q1_dot; q2_dot; q3_dot];
pd_dd=[4; 2]; %desired acceleration
J_d = diff(J, q1) * q1_dot + diff(J, q2) * q2_dot + diff(J, q3) * q3_dot;
q_dotdot=pinv(J)*(pd_dd-J_d*q_dot);

%q_dot=[0;0;0] because the robot is at rest
q_dotdot=eval(subs(q_dotdot, [q1, q2, q3, q1_dot, q2_dot, q3_dot], [pi/6, pi/6, pi/6,0,0,0]));

%we found the following acceleration command:
% q_dotdot =
% 
%     2.6524
%    -3.2466
%    -4.2175
% -> we can notice that the third link do not respect the imposed hard bound

%we can use the SNS method

qsns3_dd=-4; %we impose the violating link acceleration to its bound

psns3_dd=pd_dd-J([1 2], 3)*qsns3_dd; %we find the cartesian velocity
psns3_dd=subs(psns3_dd, [q1, q2, q3], [pi/6, pi/6, pi/6]);

%now i update the solution for the other joints
qsns12_dd=inv(J([1 2], [1 2]))*(psns3_dd-J_d*q_dot);
qsns12_dd=eval(subs(qsns12_dd, [q1, q2, q3, q1_dot, q2_dot, q3_dot], [pi/6, pi/6, pi/6,0,0,0]));

%the new command is:
qsns_dd=[qsns12_dd([1]); qsns12_dd([2]); qsns3_dd]

% qsns =
% 
%     2.7321
%    -3.4641
%    -4.0000

 %we can see that all the bounds are respected