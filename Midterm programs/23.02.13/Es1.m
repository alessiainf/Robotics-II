clc; clear all
syms l q1 al1 d1 a1 l2 q2 m q4 al2 d2 a2 l3 q3 al3 l4 d3 m4 a3 m1 m2 a1n a2n q1n q2n m3 g0 t1 t2 real
syms q1_dot q2_dot q3_dot real

% direct kinematics of 3R planar robot:
p=[l*cos(q1)+l*cos(q1+q2)+l*cos(q1+q2+q3);
    l*sin(q1)+l*sin(q1+q2)+l*sin(q1+q2+q3);
    q1+q2+q3+q4];

q=[q1;q2;q3];
q_dotdot=[q1_dot;q2_dot;q3_dot];

J=jacobian(p, q);
J=subs(J, [q1,q2,q3], [pi/4, -pi/2,pi/2]);

%the robot is at rest (c=0) and move in an horizontal plane (g=0)
%so the dynamic model is:
syms M real
u=M*q_dotdot;



%kinetic energy
w1=q1_dot;
vc1=[-l/2*sin(q1)*q1_dot; l/2*cos(q1)*q1_dot];
T1=simplify(1/2*m*norm(vc1)^2+1/2*w1'*(m*l^2/12)*w1);

w2=q1_dot+q2_dot;
vc2=[-l*sin(q1)*q1_dot-l/2*sin(q1+q2)*(q1_dot+q2_dot); 
      l*cos(q1)*q1_dot+l/2*cos(q1+q2)*(q1_dot+q2_dot)];
T2=simplify(1/2*m*norm(vc2)^2+1/2*w2'*(m*l^2/12)*w2);

w3=q1_dot+q2_dot+q3_dot;
vc3=[-l*sin(q1)*q1_dot-l*sin(q1+q2)*(q1_dot+q2_dot)-l/2*sin(q1+q2+q3)*(q1_dot+q2_dot+q3_dot); 
 l*cos(q1)*q1_dot+l*cos(q1+q2)*(q1_dot+q2_dot)+l/2*cos(q1+q2+q3)*(q1_dot+q2_dot+q3_dot)];
T3=simplify(1/2*m*norm(vc3)^2+1/2*w3'*(m*l^2/12)*w3);

T=T1+T2+T3;
T=collect(T, [q1_dot, q2_dot, q3_dot]);


q_dot=[q1_dot; q2_dot; q3_dot];
%intertia matrix
M=simplify(hessian(T,q_dot));
M=eval(subs(M, [q1,q2,q3], [pi/4, -pi/2,pi/2]))

%1)
pd=[1;0;0];


tA=M*J'*inv(J*J')*pd

%2)
tB=0;

%3)
tC=J'*inv(J*inv(M)*J')*pd