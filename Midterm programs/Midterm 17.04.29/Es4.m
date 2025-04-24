clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real

q=[q1;q2;q3];
p=[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
    sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
J=jacobian(p, q);
J=eval(subs(J, [q1, q2, q3], [pi/2, pi/3, -2*pi/3]));
vd=[1; -sqrt(3)];

%1)
H=sin(q2)^2+sin(q3)^2;

%We need to find a non singular submatrix
det([J(:, 1), J(:, 2)]); %0
det([J(:, 1), J(:, 3)]); %-1.7321
det([J(:, 2), J(:, 3)]); %-0.8660

%we choose the one with the largest determinant -> det([J1(:, 1), J1(:, 3)]) %-0.8660
%-> q1 and q3
Ja=[J(:, 1), J(:, 3)];
Jb=J(:, 2);

 
dqH=[diff(H, q1);diff(H, q2);diff(H, q3)];
dqH=eval(subs(dqH, [q1, q2, q3], [pi/2, pi/3, -2*pi/3]));

dqH_=[-(inv(Ja)*Jb)', eye(1)]*dqH;

qdb_dot=dqH_;
qda_dot=inv(Ja)*(vd-Jb*qdb_dot);

qd_dot=[qda_dot(1); qdb_dot; qda_dot(2)] %because qa is computed with q1 and q3, we need to put the correct order
% qd_dot =
% 
%    -0.4330
%     0.8660
%    -2.0000


%2)
p2=[cos(q1)+cos(q1+q2); %x2
    sin(q1)+sin(q1+q2)]; %y2
x2=cos(q1)+cos(q1+q2);
y2=sin(q1)+sin(q1+q2);

p_ta=x2^2+(y2-1.5)^2-0.75;
J_ta=simplify(jacobian(p_ta, q));
J_ta=eval(subs(J_ta, [q1, q2, q3], [pi/2, pi/3, -2*pi/3]));

%the augmented jacobian is:
JA=[J; J_ta];

%let's check for algorithmic singularities:
rank(J); %2
rank(J_ta); %1
rank(JA); %3
%since 3=2+1 -> no algorithmic singularity, 
 
% we can find the velocity command as follows:
qd_dot=inv(JA)*[vd;0]
% qd_dot =
% 
%    -0.0000
%          0
%    -2.0000



