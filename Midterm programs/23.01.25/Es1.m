clc; clear all;
syms l1 q1 al1 d1 a1 l2 q2 al2 d2 a2 l3 q3 al3 d3 a3 m1 m2 m3 g0 dc1 dc2 real
syms q1_dot q2_dot q3_dot real
syms q1_dotdot q2_dotdot q3_dotdot real
syms I1 I2 I3 real
syms rc1x rc1y rc2x rc2y rc3x rc3y real
syms a1 a2 a3 real

%A)

w1=q1_dot; 
vc1=[-dc1*sin(q1)*q1_dot; -dc1*cos(q1)*q1_dot]; 
T1=simplify(1/2*m1*norm(vc1)^2+1/2*w1'*I1*w1);

w2=q1_dot+q2_dot;
%i find the velocity of the center of mass as the derivative of the
%position of the center of mass
vc2=[-l1*sin(q1)*q1_dot;
    l1*cos(q1)*q1_dot];
T2=simplify(1/2*m2*norm(vc2)^2+1/2*w2'*I2*w2);

T=expand(T1+T2);
T=collect(T, [q1_dot, q2_dot]);

%((m1*dc1^2)/2 + (m2*l1^2)/2 + I1/2 + I2/2)*q1_dot^2 
% + I2*q1_dot*q2_dot 
% + (I2*q2_dot^2)/2


m11=2*simplify(((m1*dc1^2)/2 + (m2*l1^2)/2 + I1/2 + I2/2));

m12=simplify(I2);

m21=m12;

m22=simplify(I2);


%The inertia matrix is:
M=[m11, m12;
    m21, m22];


% %Now i find the coriolis and centrifugal terms:
q=[q1 q2];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1));
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2));

q_dot=[q1_dot; q2_dot];
c1=q_dot'*C1*q_dot;
c2=q_dot'*C2*q_dot;

c=simplify([c1;c2]); %0
 
%now i find the gravity terms
syms any real
r0c1=[dc1*cos(q1) ; any ; any];
r0c2=[l1*cos(q1) ;any ; any];

U1=simplify(-m1*[g0, 0, 0]*r0c1);
U2=simplify(-m2*[g0, 0, 0]*r0c2);
U=simplify(U1+U2);

%the gravity term is:
g_q=jacobian(U, q);
g_q=expand(simplify(transpose(g_q)));

%the viscous terms are:
syms Fv1 Fv2 real
Fv=[Fv1*q1_dot; Fv2*q2_dot];


%Now i define the dynamic terms by assuming to kno wl1 and g0:
syms a1 a2 a3 a4 a5 real
% a1= m1*dc1^2 + m2*l1^2 + I1 + I2;
% a2=I2;
% a3=dc1*m1+l1*m2;
%a4=Fv1;
%a5=Fv2;

m11=subs(m11, [m1*dc1^2 + m2*l1^2 + I1 + I2], [a1]);
m12=subs(m12, [I2], [a2]);
m21=m12;
m22=subs(m22, [ I2], [ a2]);
%The new inertia matrix is:
M=[m11, m12;
    m21, m22];

%the new gravity term is:
g_q=subs(g_q, [dc1*g0*m1*sin(q1) + g0*l1*m2*sin(q1)], [a3*sin(q1)]);

%the new viscous term is:
Fv=[a4*q1_dot; a5*q2_dot];

%Now i can build the regressor matrix as:
q_dotdot=[q1_dotdot; q2_dotdot];
u=M*q_dotdot+g_q+ Fv; %dynamic model of the robot
a=[a1;a2;a3;a4;a5]; %dynamic coefficients
Y=simplify(jacobian(u,a)); %regressor matrix



%D)
syms Itot t T real
%On the horizontal plane and without dissipative terms, the dynamic model reduces to
u=M*q_dotdot;

%%rest-to-rest trajectory while q2_dot=q2_dotdot=0:
u=subs(u, [q2_dotdot, a1, a2], [0,Itot,I2])
U1max=u(1);
U2max=u(2);

%i can use a bang-bang time profile because we aim to minimize time under a constrained force/torque:
%acceleration
q1d_dotdot_1=U1max/Itot;
q1d_dotdot_2=-U1max/Itot;

%velocity
q1d_dot_1=(U1max/Itot)*t;
q1d_dot_2=(-U1max/Itot)*(T-t);

%position
q1d_1=1/2*(U1max/Itot)*t^2;
q1d_2=1/2*(-U1max/Itot)*(T-t)^2+U1max*(T/2)^2;

%the minimum time is obtained as the area below the velocity profile
delta=U1max*(T/2)^2;
T=sqrt((4*delta*Itot)/U1max)





% Parametri simbolici/numerici
U1max = 1;    % Valore massimo del controllo (puoi modificarlo)
Itot = 1;     % Inerzia totale (modificabile)
I2=0.5;
delta = 1;    % Spostamento desiderato

T = sqrt((4*delta*Itot)/U1max);
t = linspace(0, T, 1000);

% Accelerazione bang-bang
q1d_dotdot = zeros(size(t));
q1d_dotdot(t <= T/2) = U1max / Itot;
q1d_dotdot(t > T/2) = -U1max / Itot;

% Velocità integrata
q1d_dot = cumtrapz(t, q1d_dotdot);

% Plot accelerazione
subplot(2,1,1);
plot(t, q1d_dotdot, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('q1\_dotdot');
title('Bang-bang acceleration profile');
grid on;

% Plot velocità
subplot(2,1,2);
plot(t, q1d_dot, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('q1\_dot');
title('Velocity profile from bang-bang acceleration');
grid on;


% Controllo bang-bang U1(t) e U2(t)
U1_t = Itot * q1d_dotdot;
U2_t = I2 * q1d_dotdot;

% Plot del controllo U1
figure;
subplot(2,1,1)
plot(t, U1_t, 'm', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('U1 (Torque)');
title('Bang-bang control input U1');
grid on;

% Plot del controllo U2
subplot(2,1,2)
plot(t, U2_t, 'g', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('U2 (Torque)');
title('Bang-bang control input U2');
grid on;