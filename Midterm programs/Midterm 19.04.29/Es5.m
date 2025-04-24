clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real

q=[q1; q2];

%the given inertia matrix is:
M=[a1+2*a2*cos(q2), a3+a2*cos(q2);
    a3+a2*cos(q2), a3];

%the direct kinematics is:
p=[l1*cos(q1)+l2*cos(q1+q2);
    l1*sin(q1)+l2*sin(q1+q2)];

%the distance of the robot end effect from the base is then:
p_d=simplify(norm(p));

%So we can then compute the jacobian as follow:
J=jacobian(p_d, q);

syms pd_dot real

%the velocity command can be found using the weighted jacobian based method
%with the aim to minimize the kinematic energy:
qd_dot=simplify(inv(M)*J'*inv((J*inv(M)*J'))*pd_dot);
qd_dot=eval(subs(qd_dot, [q1,q2,l1,l2,a1,a2,a3, pd_dot], [0, pi/2, 1, 1, 10, 2.5, 5/3, 0.5]))


% qd_dot =
% 
%     0.1179
%    -0.7071