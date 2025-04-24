clear all; clc;
syms q1 q2 q3 q4

ve=[0.4330; -0.75];
q=[q1, q2, q3, q4];
p=[cos(q1)+cos(q1+q2)+cos(q1+q2+q3)+cos(q1+q2+q3+q4);
    sin(q1)+sin(q1+q2)+sin(q1+q2+q3)+sin(q1+q2+q3+q4)]; %direct kinematics
Je=jacobian(p, q)
Je=eval(subs(Je, [q1, q2, q3, q4], [0, pi/6, -pi/3, -pi/3]));

vt=[-0.5; 0.8660];
p=[cos(q1)+cos(q1+q2);
   sin(q1)+sin(q1+q2)]; %direct kinematics
Jt=jacobian(p, [q1, q2])
Jt=eval(subs(Jt, [q1, q2,], [0, pi/6]));
Jt=[Jt, zeros(2,2)];


%a)
%to execute the rask velocity while minimizing q we should use the
%pseudoinverse:
disp("The velocity command required for point a) is:")
qa=pinv(Je)*ve
disp("The error for the first task is:")
ee=ve-Je*qa 
et=vt-Jt*qa %0
n=norm(et)

%b)
%to execute the rask velocity while minimizing q we should use the
%pseudoinverse:
disp("The velocity command required for point b) is:")
qb=pinv(Jt)*vt
disp("The error for the second task is:")
ee=ve-Je*qb 
et=vt-Jt*qb 
n=norm(ee)

%c)
J=[Je; Jt];
v=[ve; vt];
disp("The velocity command required for point c) is:")
qc=pinv(J)*v
disp("The error for the task is:")
e=v-J*qc
n=norm(e)


%d) 
%we should use the priority task method
Pd=eye(4)-pinv(Je)*Je;
disp("The velocity command required for point d) is:")
qd=pinv(Je)*ve+pinv((Jt*Pd))*(vt-Jt*pinv(Je)*ve)
disp("The error for the task is:")
ee=ve-Je*qd
et=vt-Jt*qd
n=norm(et)


%e) 
%we should use the priority task method
Pe=eye(4)-pinv(Jt)*Jt;
disp("The velocity command required for point e) is:")
qe=pinv(Jt)*vt+pinv((Je*Pe))*(ve-Je*pinv(Jt)*vt)
disp("The error for the task is:")
ee=ve-Je*qe
et=vt-Jt*qe
n=norm(ee)
