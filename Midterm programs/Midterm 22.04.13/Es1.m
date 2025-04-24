clear all, clc;
syms q1 q2 q3 l1 l2 l1_n l2_n

%real position of the robot
p=[l1*cos(q1)+l2*cos(q1+q2);
   l1*sin(q1)+l2*sin(q1+q2)];

%nominal position of the robot
p_n=[l1_n*cos(q1)+l2_n*cos(q1+q2);
   l1_n*sin(q1)+l2_n*sin(q1+q2)];

%delta_p=p-p_n=[delta_l1*cos(q1)+delta_l2*cos(q1+q2); = FI_q*delta_l
%               delta_l1*sin(q1)+delta_l2*sin(q1+q2)];

%where:
FI_q=[cos(q1), cos(q1+q2);
       sin(q1), sin(q1+q2)];
%delta_l=[delta_l1; delta_l2];

%in order to find the values that we are looking for we need to compute the
%following values:
%l1=l1_n+delta_l1;
%l2=l2_n+delta_l2;

%we first need to find the value of delta_l

%nominal positions of the robot using data:
p_n=subs(p_n, [l1_n, l2_n], [1, 1]);
p_na=eval(subs(p_n, [q1, q2], [0, 0]));
p_nb=eval(subs(p_n, [q1, q2], [pi/2, 0]));
p_nc=eval(subs(p_n, [q1, q2], [pi/4, -pi/4]));
p_nd=eval(subs(p_n, [q1, q2], [0, pi/4]));

%matrix of the measures
delta_p=[[2;0]-p_na; [0;2]-p_nb; [1.6925; 0.7425]-p_nc; [1.7218;0.6718]- p_nd];

FI_a=eval(subs(FI_q, [q1, q2], [0, 0]));
FI_b=eval(subs(FI_q, [q1, q2], [pi/2, 0]));
FI_c=eval(subs(FI_q, [q1, q2], [pi/4, -pi/4]));
FI_d=eval(subs(FI_q, [q1, q2], [0, pi/4]));

FI_q=[FI_a; FI_b; FI_c; FI_d];

delta_l=pinv(FI_q)*delta_p; % ==[delta_l1; delta_l2]

disp("the values we are looking for are:")
l1_n=1; l2_n=1; 
delta_l1=delta_l(1); delta_l2=delta_l(2);
l1=l1_n+delta_l1
l2=l2_n+delta_l2



