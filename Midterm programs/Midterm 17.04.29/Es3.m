clc; clear all
clc; clear all;
syms l q1 al1 d1 a1 L2 q2 q4 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 a1n a2n q1n q2n m3 g0 t1 t2 real
syms q1_dot q2_dot q3_dot real
syms rc1x rc2x rc2y rc3x rc3y real% center of mass not on the link axis
syms Ic1xx Ic1yy Ic1zz Ic2xx Ic2yy Ic2zz Ic3xx Ic3yy Ic3zz real

% direct kinematics of 4R planar robot:
p=[l*cos(q1)+l*cos(q1+q2)+l*cos(q1+q2+q3)+l*cos(q1+q2+q3+q4);
    l*sin(q1)+l*sin(q1+q2)+l*sin(q1+q2+q3)+l*sin(q1+q2+q3+q4)];

p=subs(p, l, 0.5);
q=[q1;q2;q3;q4];

%initial velocity
v_in=[0;10];


%the velocity command is:
J=jacobian(p, q);
%initial configuration
J=subs(J, [q1, q2, q3, q4], [0, 0, 0, 0]);

%minimum norm solution
q_dot=eval(pinv(J)*v_in);

%we can see that the third joint is over the limit. 
% q_dot =
% 
%     2.6667
%     2.0000
%     1.3333
%     0.6667

%we can use the SNS method
q_dot_sns3=1; %saturate the violating joint

%recompute the cartesian velocity
p_dot_sns=v_in-J(1:2, 3)*q_dot_sns3;

%update the solution for all the joints
J_j=[J(1:2, 1), J(1:2, 2), J(1:2, 4)];
q_dot_sns124=eval(pinv(J_j)*p_dot_sns);

q_dot_sns=([q_dot_sns124(1); q_dot_sns124(2); q_dot_sns3; q_dot_sns124(3)]);

%the second joint violates the hard bound
% q_dot_sns =
% 
%     2.7692
%     2.0769
%     1.0000
%     0.6923
%we can use again the SNS method


q_dot_sns2=2; %saturate the violating joint

%recompute the cartesian velocity
p_2_dot_sns=v_in-J(1:2, 2)*q_dot_sns2-J(1:2, 3)*q_dot_sns3;

%update the solution for all the joints
J_2_j=[J(1:2, 1), J(1:2, 4)];
q_dot_sns14=eval(pinv(J_2_j)*p_2_dot_sns);

q_dot_sns=([q_dot_sns14(1); q_dot_sns2; q_dot_sns3; q_dot_sns14(2)]);

%The bound are respected and the final solution is:
% q_dot_sns =
% 
%     2.8235
%     2.0000
%     1.0000
%     0.7059
