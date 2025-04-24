clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q4 d4 q3 al3 d3 a3 m1 m2 m3 m4 g0 dc1 dc2 l real
syms q1_dot q2_dot q3_dot q4_dot real
syms q1_dotdot q2_dotdot q3_dotdot q4_dotdot real


%direct kinematics
p=[q1+l*cos(q3);
q2+l*sin(q3)];
q=[q1; q2; q3];
J=jacobian(p, q);
rank(J); %2 -> we can use jacobian based method
J_p=pinv(J);

%1)
%m, m/s
p_d=[-1; 1];
q_dot=J_p*p_d;
q_dot=eval(subs(q_dot, [l, q3],[0.5, pi/6]));
% q_dot =
% 
%    -0.8634
%     0.7634
%     0.5464


%cm, cm/s
p_d=[-100; 100];
q_dot=J_p*p_d;
q_dot=eval(subs(q_dot, [l, q3],[50, pi/6]));
q_dot=[q_dot(1)*0.01; q_dot(2)*0.01; q_dot(3)];
% q_dot =
% 
%    -0.3173
%    -0.1825
%     2.7310

%this happen because det(J*J')=1+l^2 -> the pseudoinverse computation depends on the units chosen
deter=simplify(det(J*J'));

%2)
syms w real
W=[1,0,0;
    0,1,0;
    0,0,w];
Jw_p=inv(W)*J'*inv((J*inv(W)*J'));

deter=simplify((det(J*inv(W)*J'))); %(l^2 + w)/w
% -> for w=l^2 we can remove the problematic term

%m, m/s
p_d=[-1; 1];
q_dot=Jw_p*p_d;
q_dot=eval(subs(q_dot, [l, q3, w],[0.5, pi/6, 0.5^2]));
% q_dot =
% 
%    -0.6585
%     0.4085
%     1.3660

%cm, cm/s
p_d=[-100; 100];
q_dot=Jw_p*p_d;
q_dot=eval(subs(q_dot, [l, q3, w],[50, pi/6, 50^2]));
q_dot=[q_dot(1)*0.01; q_dot(2)*0.01; q_dot(3)];
% q_dot =
% 
%    -0.6585
%     0.4085
%     1.3660

%3)
%The use of very large weights slow down the movements of the joints. 
%So if we use a very large weight on the third joint, only the first two
%joints will move -> w=100000

%m, m/s
p_d=[-1; 1];
q_dot=Jw_p*p_d;
q_dot=eval(subs(q_dot, [l, q3, w],[0.5, pi/6, 100000]));
% q_dot =
% 
%    -1.0000
%     1.0000
%     0.0000