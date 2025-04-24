clc; clear all;
syms l1 L q1 al1 d1 a1 l2 q2 al2 d2 q4 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot q4_dot real
syms q1_dotdot q2_dotdot q3_dotdot q4_dotdto real 
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real
syms a1 a2 a3 real

q=[q1 q2 q3 q4];

%task with priority
p1=[L*cos(q1)+L*cos(q1+q2)+(L/4)*cos(q1+q2+q3)+(L/4)*cos(q1+q2+q3+q4);
   L*sin(q1)+L*sin(q1+q2)+(L/4)*sin(q1+q2+q3)+(L/4)*sin(q1+q2+q3+q4)];
J1=jacobian(p1, q);

%task with second priority
p2=[L*cos(q1)+L*cos(q1+q2);
   L*sin(q1)+L*sin(q1+q2)];
J2=jacobian(p2, q);


%the extended jacobian is:
Je=[J1;J2];
det=simplify(det(Je))  %(L^4*sin(q2)*sin(q4))/16

%let's notice that:
% J11=subs(J1, [q2,q4] [0,0])
% rank(J11) %the rank do not change in the algorithmic singularity
%-> the first task can be executed anyway

% J22=subs(J2, [q2,q4] [0,0])
% rank(J22) %the rank do change in the algorithmic singularity
%in particular it change for sin(q2)=0
%the second task cannot be executed when sin(q2)=0 (??)

%-> a true algorithmic singularity occurs when the extended jacobian has a
%singularity, so when sin(q2)!=0 and sin(q4)=0
 
%Under these condition (sin(q4)=0) we can still execute the first priority
%task and the second:
syms p1_dot p2_dot real

PA2=(eye(4)-pinv(J1)*J1);
qtp_dot=pinv(J1)*p1_dot+pinv(J2*PA2)*(p2_dot-J2*pinv(J1)*p1_dot);

%when also sin(q2)=0 the secondo priority task is not realized





