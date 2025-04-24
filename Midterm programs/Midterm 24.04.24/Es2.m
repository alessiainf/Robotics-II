clc; clear all;
syms g0 al1 d1 q1 al2 d2 a1 a2 q2 l1 l2 m2 real
syms q1_dot q2_dot q1_dotdot q2_dotdot real
syms Ic2xx Ic2yy Ic2zz rc1x rc1y rc1z rc2x m1 Ic1xx Ic1xy Ic1xz Ic1yx Ic1yy Ic1yz Ic1zx Ic1zy Ic1zz real

Ic1=[Ic1xx, Ic1xy, Ic1xz;
    Ic1yx, Ic1yy, Ic1yz;
    Ic1zx, Ic1zy, Ic1zz];
rc1=[rc1x; rc1y; rc1z];

Ic2=[Ic2xx, 0, 0;
    0, Ic2yy, 0;
    0, 0, Ic2zz];
rc2=[rc2x; 0; 0];

%i)
H_1 = [
    cos(q1), -sin(q1)*cos(al1), sin(q1)*sin(al1), a1*cos(q1);
    sin(q1), cos(q1)*cos(al1), -cos(q1)*sin(al1), a1*sin(q1);
    0, sin(al1), cos(al1), d1;
    0, 0, 0, 1
];

H_2 = [
    cos(q2), -sin(q2)*cos(al2), sin(q2)*sin(al2), a2*cos(q2);
    sin(q2), cos(q2)*cos(al2), -cos(q2)*sin(al2), a2*sin(q2);
    0, sin(al2), cos(al2), d2;
    0, 0, 0, 1
];

H1_ = subs(H_1, [al1, a1, d1, q1], [-pi/2, 0, l1, q1]);
R0_1=H1_(1:3, 1:3);
% %disp('punto 0p1:');
r0_01=H1_(1:3, end);


H2_ = subs(H_2, [al2, a2, d2, q2], [0, l2, 0, q2]);
R1_2=H2_(1:3, 1:3);
% %disp('punto 1p2:');
r1_12=H2_(1:3, end);

p=H1_*H2_*[0;0;0;1];
p=p(1:3);

w1=transpose(R0_1)*(0+q1_dot*[0;0;1]);
v1=transpose(R0_1)*(0+cross((R0_1*w1), r0_01));
vc1=v1+cross(w1,rc1);
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*Ic1*w1);

w2=transpose(R1_2)*(w1+q2_dot*[0;0;1]);
v2=transpose(R1_2)*(v1+cross((R1_2*w2), r1_12));
vc2=v2+cross(w2,rc2);
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*Ic2*w2);


T=expand(T1+T2);
T=collect(T, [q1_dot, q2_dot]);

q_dot=[q1_dot; q2_dot];
M=simplify(hessian(T, q_dot));
M=collect(M, cos(q2));

%dynamic coefficients
% a1=m2*l2^2 + 2*m2*l2*rc2x + m2*rc2x^2 - Ic2xx + Ic2yy
% a2=m1*rc1x^2 + m1*rc1z^2 + Ic2xx + Ic1yy
% a3=m2*l2^2 + 2*m2*l2*rc2x + m2*rc2x^2 + Ic2zz
syms a1 a2 a3 real
m11=a1*cos(q2)^2 + a2;
m12=0;
m21=m12;
m22=a3;

M=[m11, m12; 
    m21, m22]

% %Now i find the coriolis and centrifugal terms:
q=[q1 q2];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));

q_dot=[q1_dot; q2_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;

c=simplify([c1;c2]) 

S1=q_dot'*C1;
S2=q_dot'*C2;
S=[S1;S2]
 
%now i find the gravity terms
syms any real
r0c1=H1_*[rc1;1];
r0c2=H1_*H2_*[rc2;1];

U1=simplify(-m1*[-g0, 0, 0]*r0c1([1 2 3]));
U2=simplify(-m2*[-g0, 0, 0]*r0c2([1 2 3]));
U=simplify(U1+U2);

%the gravity term is:
g_q=jacobian(U, q);
g_q=expand(simplify(transpose(g_q)));
g_q=collect(g_q, [cos(q1),sin(q1), sin(q2)]);

%dynamic terms
% a4=-m1*rc1z
% a5=-*l2*m2-m2*rc2x
%a6=-*m1*rc1x
syms a4 a5 a6 real
g_q=g0*[a4*cos(q1)+a6*sin(q1)+a5*cos(q2)*sin(q1);
    a5*cos(q1)*sin(q2)]

%ii)
q_dotdot=[q1_dotdot; q2_dotdot];
u=M*q_dotdot+c+g_q; %dynamic model of the robot
a=[a1;a2;a3;a4;a5;a6]; %dynamic coefficients
Y=simplify(jacobian(u,a)) %regressor matrix

%iii)
syms t real
u1=subs(u, [q1, q2, q1_dot, q2_dot, q1_dotdot, q2_dotdot], [2*t, pi/4, 2, 0, 0, 0])

%iv)
u2=subs(u, [q1_dot, q2_dot, q1_dotdot, q2_dotdot], [ 0, 0, 0, 0])
%i need to find values of q1,q2 s.t. u2=0:

% cos(q1)=0 -> q1=+- pi/2
% -> | +- (a6* + a5*cos(q2)) |=0  -> cos(q2)= +- (a6/a5) -> q2=+-arccos(a6/a5)
%    |           0           |=0

%sin(q2)=0 -> q2=0 | pi
% -> | a4*cos(q1) + a6*sin(q1) +-a5*sin(q1)) |=0  -> 
%    |           0                           |=0

% -> | a4*cos(q1) + (a6 +- a5)*sin(q1) ) |=0  -> cos(q1)/sin(q1)= -(a6 +-a5) /a4
%    |           0                       |=0

% -> q1=arccotan (-(a6 +-a5) /a4)


%v)
% g_q =
% g0*(a4*cos(q1) + a6*sin(q1) + a5*cos(q2)*sin(q1))
%                             a5*g0*cos(q1)*sin(q2)

% a4=-m1*rc1z =0 <-> rc1z =0
% a5=-*l2*m2-m2*rc2x=0 <-> rc2x=-l2
%a6=-*m1*rc1x=0 <-> rc1x=0
 
%vi)
%i need to find values of q1,q2 s.t. u2=0:

%robot at rest c=0 and g(q)=0 for previous conditions
%Jyz'*Fyz=M*q_dotdot
%we're interested in y and z direction
%qyz_dotdot=Jyz'*p_dotdot

%-> Jyz'*Fyz=M*Jyz'*p_dotdot

%-> p_dotdot=Jyz*inv(M)*Jyz'*Fyz

% Mp_inv=J*inv(M)*Jyz'
% 
% p_dotdot=Mp_inv*Fyz;
syms Fy Fz real
J=jacobian(p, q);
Jyz=J([2 3], [1 2]);
Fyz=[Fy; Fz];

Mp_inv=Jyz*inv(M)*Jyz' %-> not a diagonal matrix -> the acceleration will NOT have the same direction as the tip
p_dotdot=Mp_inv*Fyz
