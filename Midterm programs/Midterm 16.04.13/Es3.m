clc; clear all
clc; clear all;
syms l q1 al1 d1 a1 l2 q2 q4 al2 d2 a2 l3 q3 al3 l4 d3 m4 a3 m1 m2 a1n a2n q1n q2n m3 g0 t1 t2 real
syms q1_dot q2_dot q3_dot real

% direct kinematics of 4R planar robot:
p=[l*cos(q1)+l*cos(q1+q2)+l*cos(q1+q2+q3)+l*cos(q1+q2+q3+q4);
    l*sin(q1)+l*sin(q1+q2)+l*sin(q1+q2+q3)+l*sin(q1+q2+q3+q4);
    q1+q2+q3+q4];

q=[q1;q2;q3;q4];

J=jacobian(p, q);
rank(J); % 3

detJL=simplify(det(J*J'));
%-l^4*(cos(q2 + 2*q3) + cos(2*q2 + q3) + cos(2*q2) + cos(2*q3) + cos(2*q2 + 2*q3) - cos(q2) - cos(q3) - 3)
%detJL=0 if q2=0 or pi, q3=0 or pi

%singular configurations:
%J=subs(J, [q2,q3], [0,0]) %rank=2
%J=subs(J, [q2,q3], [pi,pi]) %rank=2

vd=[1;0];
fid_dot=0.5;
pd_dot=[vd;fid_dot];
J=subs(J, [l, q1,q2,q3,q4], [0.5, 0,0,pi/2,0]);
rank(J); %rank=3 maximum

%range function:
H=(1/(2*4))*((q1/4)^2+(q2/4)^2+(q3/4)^2+(q4/4)^2);
%gradient of the range function:
dH=[diff(H, q1);diff(H, q2);diff(H, q3);diff(H, q4)];
dH=eval(subs(dH, [q1,q2,q3,q4], [0,0,pi/2,0]));

%velocity command using the projected gradient method:
q_dot=eval(pinv(J)*pd_dot+(eye(4)-pinv(J)*J)*(-dH))







