clc; close all; clear;
syms q1 q2 q3 real

p=[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
    sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
q=[q1 q2 q3];
J=jacobian(p, q);
J=eval(subs(J, [q1, q2, q3], [2*pi/5, pi/2, -pi/4]));
J_p=pinv(J);
rd=[-3;0]; %vx

%The H function is the Hrange:
H1=(q1/(pi/2+pi/2))^2;
H2=((q2-((2*pi/3-0)/2))/(2*pi/3))^2;
H3=(q3/(pi/4+pi/4))^2;
H=1/6 *(H1+H2+H3);

dH=gradient(H, q);
dH=eval(subs(dH, [q1, q2, q3], [2*pi/5, pi/2, -pi/4]));

%the velocity command is found with the jacobian based method:
qd_dot=J_p*rd-(eye(3)-J_p*J)*dH %notice that in Hrange, we use -dh and not dh

% qd_dot =
%     2.0638
%    -1.9261
%     0.9786
%we can see that the velocity exceed in the first component
%so we need to scale it, in particular we need to slow it down.

qr=J_p*rd;

% qr =
%     2.1076
%    -1.9261
%     0.8730

qn=-(eye(3)-J_p*J)*dH;

% qn =
%    -0.0437
%    -0.0000
%     0.1056

%the null space part do not change, we need to change the first part that
%regards the minimum norm solution
%k=(qmax-q1n)/q1r
k=(2-(-0.0437))/2.1076

%the new solution is:
qd_dot=k*qr+qn;

% qd_dot =
%     1.9999
%    -1.8677
%     0.9521


%note: A direct application of the SNS method to recover feasibility would not be correct, 
% since the solution qd_dot contains a null-space term that does not scale with the task 
% velocity r_dot .