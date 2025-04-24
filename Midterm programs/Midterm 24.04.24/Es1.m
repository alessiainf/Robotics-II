% Robotics 2 Midterm
% 24 April 2024
% Ex #1
%
% Projected Gradient (PG) method (with obstacle avoidance) for a 3R planar robot
clear all, clc;
syms q1 q2 q3 real

%i)
q=[q1; q2; q3];
pe_simb=[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);sin(q1)+sin(q1+q2)+sin(q1+q2+q3)]; %Direct kinematics
Je_simb=jacobian(pe_simb,q);
q0=[0;pi/2;-pi/2];
pe=subs(pe_simb,q,q0);
Je=subs(Je_simb,q,q0);
Je_pinv=pinv(Je);
ve=[0; 1];
%The function H is the one that maximizes the distance between the robot
%and the obstacle so it can be written as follow:
C=[0;2];r=0.5;
pm_simb=[cos(q1)+cos(q1+q2);sin(q1)+sin(q1+q2)]; %Direct kinematics of the nearest point
pm=subs(pm_simb,q,q0); % (1,1)
Jm_simb=jacobian(pm_simb,q);
Jm=subs(Jm_simb,q,q0);
dH=0.5*Jm'*(pm-C)/norm(pm-C);

disp("The value of the velocity command is the following:");
dq_PG=Je_pinv*ve+(dH-Je_pinv*Je*dH); %joint velocity
disp(eval(dq_PG));

disp("In the given configuration the robot nearest point is in (1,1) of the joint 2")
disp("So we can find the velocity of vm as the velocity of joint 2:")
vm=eval(Jm*dq_PG);  %task velocity
disp(vm);

%ii)
% first case
Pe=eye(3)-Je_pinv*Je;
dq_TP_A=Je_pinv*ve+pinv(Jm*Pe)*(vm-Jm*Je_pinv*ve);
disp("The velocity with the task priority method is the following:")
disp(eval(dq_TP_A));
disp("The joint velocities of the two methods have the same directions.")

% second case
alfa=1;
vm_B=alfa*(1-(r/norm(pm-C)))*(pm-C);
dq_TP_B=Je_pinv*ve+pinv(Jm*Pe)*(vm_B-Jm*Je_pinv*ve);
disp("The velocity with the task priority method is the following:")
disp(eval(dq_TP_B));
disp("The joint velocities of the two methods do not have the same directions.")