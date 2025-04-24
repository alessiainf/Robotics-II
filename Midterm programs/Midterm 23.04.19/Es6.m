clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real
syms a1 a2 a3 real

%i first find the inertia matrix
w1=0; %prismatic joint
vc1=q1_dot;
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*I1*w1);
w2=q2_dot;
vc2=[q1_dot; dc2*q2_dot];
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*I2*w2);
T=expand(T1+T2);
T=collect(T, [q1_dot, q2_dot]);
q_dot=[q1_dot; q2_dot];
%intertia matrix
M=simplify(hessian(T,q_dot));

%Now i find the coriolis and centrifugal terms:
q=[q1 q2];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));
q_dot=[q1_dot; q2_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;
c=simplify([c1;c2]);

%Now i find the gravity terms
syms any real
r0c1=[any; q1 ; any];
r0c2=[any ;q1; any];
U1=simplify(-m1*[0, -g0, 0]*r0c1);
U2=simplify(-m2*[0, -g0, 0]*r0c2);
U=simplify(U1+U2);
%the gravity term is:
g_q=jacobian(U, q);
g_q=expand(simplify(transpose(g_q)));

%The dynamic model of the robot is u:
q_dotdot=[q1_dotdot; q2_dotdot];
u=M*q_dotdot+c+g_q; 


%the rest-to-rest cubic trajectory is qd. The acceleration is:
syms t T delta_q1 delta_q2 real
delta_q=[delta_q1; delta_q2];
qd_dotdot=(delta_q/T^2)*(6-12*(t/T));

%Now i substitute qd_dotdot in u and consider that the maximum torque would be at the
%start (t=0) or the end (t=T) of the trajectory:
ud=(subs(u, q_dotdot, qd_dotdot));
ud=subs(ud, t, 0);

% ud =
% g0*m1 + g0*m2 + (6*delta_q1*(m1 + m2))/T^2   -> i can extraxt T1
%           (6*delta_q2*(m2*dc2^2 + I2))/T^2   -> i can extract T2


%So now i can find T from u
syms u1_max u2_max real
% g0*m1 + g0*m2 + (6*delta_q1*(m1 + m2))/T^2 -> 
T1=sqrt( (g0*m1 + g0*m2 + (6*delta_q1*(m1 + m2))) / u1_max);
% (6*delta_q2*(m2*dc2^2 + I2))/T^2 ->
T2=sqrt( (6*delta_q2*(m2*dc2^2 + I2)) / u2_max);

%minimun feasible time
T=min(T1, T2)


