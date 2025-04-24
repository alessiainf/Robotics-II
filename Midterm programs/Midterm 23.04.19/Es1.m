clc; clear; close all;
% Definizione delle variabili simboliche
syms q1 q2 q3 q4 l px py real % Variabili generalizzate (angoli e distanze)
syms q1_dot q2_dot q3_dot real % Velocità articolari
syms q1_ddot q2_ddot q3_ddot real    % Accelerazioni articolari

%Direct kinematics of the robot
p = [l*cos(q1)+l*cos(q1+q2)+l*cos(q1+q2+q3);
    l*sin(q1)+l*sin(q1+q2)+l*sin(q1+q2+q3)];

q = [q1, q2, q3]; 
q_dot = [q1_dot; q2_dot; q3_dot];


%analytic jacobian
J=simplify(jacobian(p, q));

%Cartesian acceleration is: p_dotdot=J*q_dotdot+n where: n=J_dot*q_dot
J_dot = diff(J, q1) * q1_dot + diff(J, q2) * q2_dot + diff(J, q3) * q3_dot;
n=J_dot*q_dot;
%inverting this equation we can find the value of the joint acceleration
pe_dotdot=[2; 1]; %desired cartesian acceleration
q_dotdot=pinv(J)*(pe_dotdot-n);
q_dotdot=eval(subs(q_dotdot, [q1, q2, q3, q1_dot, q2_dot, q3_dot, l], [0, 0, pi/2, 0.8, 0, -0.8, 0.5]))


%Now we need to verify if this command is feasible by considering both the
%acceleration and velocity bounds.
Tc=0.1;

%in the time interval the acceleration is kept costant. The velocity is:

%q_dot_next=q_dot+q_dotdot*Tc -> (q_dotdot=q_dot_next-q_dot)/Tc

% ->where we impose q_dot_next=v_max to see the feasible maximum joint acceleration

A1max=10; V1max=1.5; Av1max=(V1max+0.8)/Tc; Av1min=(V1max-0.8)/Tc;
A2max=10; V2max=1.5; Av2max=(V2max+0)/Tc;   Av2min=(V2max-0)/Tc;
A3max=10; V3max=1;   Av3max=(V3max-0.8)/Tc; Av3min=(V3max+0.8)/Tc;

Q1min=max(-A1max, -Av1max); Q1max=min(A1max, Av1min);
Q2min=max(-A2max, -Av2max); Q2max=min(A2max, Av2min);
Q3min=max(-A3max, -Av3max); Q3max=min(A3max, Av3min);

Qmin=[Q1min; Q2min; Q3min]
Qmax=[Q1max; Q2max; Q3max]

%we notice that the third joint exceeds the lower limit (-5.4 < -2):
% q_dotdot =
%     1.8800
%    -1.7600
%    -5.4000
% 
% Qmin =
%   -10.0000
%   -10.0000
%    -2.0000
% 
% Qmax =
%     7.0000
%    10.0000
%    10.0000


%we can then apply the SNS algorithm at the acceleration level to redifine
%the acceleration command

%i find the joint that do not respect the bound and update it to the bound:
ddqsns3=-2; 

%i recompute the cartesian acceleration
ddpsns3=pe_dotdot-(J([1 2], [3])*ddqsns3);
ddpsns3=subs(ddpsns3, [l, q1, q2, q3], [0.5, 0, 0, pi/2]);

%i find the value of the other joint's accelerations 
ddqsns12=inv(J([1 2], [1 2]))*(ddpsns3-n);
ddqsns12=eval(subs(ddqsns12, [l, q1, q2, q3, q1_dot, q2_dot, q3_dot], [0.5, 0, 0, pi/2,  0.8, 0, -0.8]));

%the new velocity command is:
ddqsns=[ddqsns12(1); ddqsns12(2); ddqsns3]
%limits are respected so we can stop here