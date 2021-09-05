clear all
close all
clc

%% Parameters
v_ref = 0.1485;             % [V]        Input Voltage 
g = 9.81;                   % [m/sec^2]  Gravity acceleration 
m = 0.03;                   % [kg]       Wheel weight
R = 0.04;                   % [m]        Wheel radius
J_w = (m*R^2)/2;            % [kg∙m^2]   Wheel inertia moment
M = 0.6;                    % [kg]       Body weight
W = 0.14;                   % [m]        Body width
D = 0.04;                   % [m]        Body depth
H = 0.144;                  % [m]        Body height
L = H/2;                    % [m]        Distance of the center of mass from the wheel axle
J_psi = (M*L^2)/3;          % [kg∙m^2]   Body pitch inertia moment
J_phi = M*((W^2)+(D^2))/12; % [kg∙m^2]   Body yaw inertia moment
J_m = 1e-5;                 % [kg∙m^2]   DC motor inertia moment
R_m = 6.69;                 % [Ω]        DC motor resistance
K_b = 0.468;                % [Vsec/rad] DC motor back EMF constant
K_t = 0.317;                % [N∙m/A]    DC motor torque constant
n = 1;                      %            Gear ratio
f_m = 0.0022;               %            Friction coefficient between body and DC motor
f_w = 0;                    %            Friction coefficient between wheel and floor

%% Symbolic Equations

syms psi psi_dot th_l th_r th_ml th_mr thl_dot thr_dot

th_dot=(thl_dot+thr_dot)/2;

% Motion Equations
th = (th_l+th_r)/2;         phi = (R/W)*(th_r-th_l);
x_m = R*th*cos(phi);        y_m = R*th*sin(phi);        z_m = R;
x_l = x_m-(W/2)*sin(phi);   y_l = y_m+(W/2)*cos(phi);   z_l = z_m;
x_r = x_m+(W/2)*sin(phi);   y_r = y_m-(W/2)*cos(phi);   z_r = z_m;
x_b = x_m+L*sin(psi)*cos(phi);  y_b = y_m+L*sin(psi)*sin(phi);  z_b = z_m+L*cos(psi);

% Energy Equations
% 1. Translation Kinetic Energy
T_11 = 0.5*m*(((diff(x_l))^2)+((diff(y_l))^2)+((diff(z_l,psi))^2));
T_12 = 0.5*m*(((diff(x_r))^2)+((diff(y_r))^2)+((diff(z_r,psi))^2));
T_13 = 0.5*M*(((diff(x_b))^2)+((diff(y_b))^2)+((diff(z_b))^2));
T_1 = T_11 + T_12 + T_13;
T_1 = vpa(T_1,4);

% 2. Rotational Kinetic Energy
% T_21 = 0.5*J_w*((diff(th_l)^2)+(diff(th_r)^2));
% T_22 = 0.5*(J_psi*((diff(psi)^2))+J_phi*((diff(phi)^2)));
% T_23 = 0.5*(n^2)*J_m*(((diff(th_l)-diff(psi))^2)+(((diff(th_r)-diff(psi))^2))); % rotation kinetic energy of an armature in left and right DC motor.
% T_2 = T_21 + T_22 + T_23;
% T_2 = double(T_2)

% 2. Rotational Kinetic Energy
T_21 = 0.5*J_w*((thl_dot^2)+(thr_dot^2));
T_22 = 0.5*(J_psi*(psi_dot^2)+J_phi*(psi_dot^2));
T_23 = 0.5*(n^2)*J_m*(((thl_dot-psi_dot)^2)+((thr_dot-psi_dot)^2)); % rotation kinetic energy of an armature in left and right DC motor.
T_2 = T_21 + T_22 + T_23;
T_2 = simplify(T_2);
T_2 = vpa(T_2,4);

% 3. Potential Energy
U = m*g*z_l + m*g*z_r + M*g*z_b;
U = vpa(U,5);

% 4. Lagrangian
Lag = T_1+T_2-U;
Lag = simplify(Lag);
Lag = vpa(Lag,4)

% Equation of forces are mentioned in the report
%% Linearized System

alpha = (n*K_t)/R_m;
beta = (n*K_t*K_b)/R_m + f_m;

% Equation Matrices
E = [(2*m + M)*R^2 + 2*J_w + 2*(n^2)*J_m,    M*L*R - 2*(n^2)*J_m
     M*L*R-2*n^2*J_m,                        M*L^2 + J_psi + 2*n^2*J_m];
F = [beta+f_m -beta
     -2*beta  2*beta];
G = [0 0
     0 -M*g*L];
H = [alpha/2 alpha/2
     -alpha -alpha];
I = beta+(W/R)*f_w;
J = (R/W)*alpha;
K = 0.5*m*W^2 + J_phi + ((W^2)/(2*R^2))*(J_w+J_m*n^2);

%% System Matrices

A1_32=-((M*g*L*E(1,2))/det(E)); 
A1_42= (M*g*L*E(1,1))/det(E); 
A1_33=-(((beta + f_w)*E(2,2)+2*beta*E(1,2)))/det(E); 
A1_43= (((beta + f_w)*E(1,2)+2*beta*E(1,1)))/det(E); 
A1_34= (beta*(E(2,2))+2*E(1,2))/det(E); 
A1_44= -(beta*(E(1,2))+2*E(1,1))/det(E); 
B1_3= ((alpha*(E(2,2)/2))+(E(1,2)))/det(E); 
B1_4= -alpha*(E(1,2)/2+ E(1,1))/det(E);


A1 = [0 0 1 0
      0 0 0 1
      0 A1_32 A1_33 A1_34
      0 A1_42 A1_43 A1_44]

B1 = [0
      0
      B1_3
      B1_4]

% A2 = [0 1;0 -I/K];
% B2 = [0 0;-J/K J/K];
C1=[0 1 0 0]; 
C=[1 1 0 0];

D=0;

%% Original System

my_sys = ss(A1,B1,C,D);
tran_sys = tf(my_sys)
eig_val = eig(my_sys);
cpoly_orig = poly(A1);

x10 = [0;(0/180)*pi;0;0];

x0 = [0;(15/180)*pi;0;0];

% Controllable
P = ctrb(my_sys);
P_det = det(P);
P_rank = rank(P);
if (rank(P) == size(A1,1)) % Logic to assess controllability
    disp('System is controllable.');
else
    disp('System is NOT controllable.');
end

% Observable
Q = obsv(my_sys);
Q_rank = rank(Q);
Q_det = det(Q);
if (rank(Q) == size(A1,1)) % Logic to assess observability
    disp('System is observable.');
else
    disp('System is NOT observable.');
end

% Plots
t=0:0.05:5;
[y,t,x]=initial(my_sys,x0,t);


%% Controllabilty Analysis

my_pole=-1+0.75i;
d_pole_c = [my_pole,5*real(my_pole),conj(my_pole),10*real(my_pole)]; %using 1.5% overshoot and 4 sec settling time
K_acker = acker(A1,B1,d_pole_c) % Controller Gain

ctrl_sys= ss((A1- B1*K_acker),B1,C,D);

A_cont = A1- B1*K_acker;
eig_A_cont = eig(A_cont)

[y_c,t_c,x_c]=initial(ctrl_sys,x0,t);

K_acker_nl = [40 -2 0 -2.1]

%% Observability Analysis

d_pole_o = [my_pole,5*real(my_pole),conj(my_pole),10*real(my_pole)]*10;
ob_acker = acker(A1',C',d_pole_o) % Observable Gain
ob_acker = ob_acker';
A_ob = A1-ob_acker*C;

obsv_sys = ss(A_ob,B1,C,D);

[y_o,t_o,x_o]=initial(obsv_sys,x0,t);

%% Plots

% Original System Plots
figure (1)
subplot(2,2,1)
plot(t,x(:,2)*(180/pi))
title('Initial System Output reponse for Psi')
ylabel('Angular Displacement (deg)')
xlabel('Time (s)')

subplot(2,2,2)
plot(t,x*(180/pi))
title('All States Initial System Response')
ylabel('Angular Displacement and Velocities')
xlabel('Time (s)')
legend('theta','psi','psi-dot','theta-dot')

subplot(2,2,3)
step(my_sys)
title('Step Response: No Controller Designed')

subplot(2,2,4)
pzmap(my_sys)
title('Initial System Pole-Zero map')

% Controllability Plots
figure (2)
subplot(2,2,1)
plot(t_c,x_c(:,2)*(180/pi))
title('Controllable System Output reponse for Psi')
ylabel('Angular Displacement (deg)')
xlabel('Time (s)')

subplot(2,2,2)
plot(t_c,x_c*(180/pi))
title('All States Controllable System Response')
ylabel('Angular Displacement and Velocities')
xlabel('Time (s)')
legend('theta','psi','psi-dot','theta-dot')

subplot(2,2,3)
step(ctrl_sys)
title('Step Response: Controller Designed')

subplot(2,2,4)
pzmap(ctrl_sys)
title('Controller System Pole-Zero map')

% Observability Plots
figure (3)
subplot(2,2,1)
plot(t_o,x_o(:,2)*(180/pi))
title('Observable System Output reponse for Psi')
ylabel('Angular Displacement (deg)')
xlabel('Time (s)')

subplot(2,2,2)
plot(t_o,x_o*(180/pi))
title('All States Observable System Response')
ylabel('Angular Displacement and Velocities')
xlabel('Time (s)')
legend('theta','psi','psi-dot','theta-dot')

subplot(2,2,3)
step(obsv_sys)
title('Step Response: Observable Designed')

subplot(2,2,4)
pzmap(obsv_sys)
title('Observable System Pole-Zero map')

% Psi State plots
figure (4)
subplot(2,2,[1,2])
plot(t,x(:,2)*(180/pi))
title('State Response for Original System')
ylabel('Psi (deg)')
xlabel('Time (s)')

subplot(2,2,3)
plot(t_c,x_c(:,2)*(180/pi))
title('State Response for Controllable System')
ylabel('Psi (deg)')
xlabel('Time (s)')

subplot(2,2,4)
plot(t_o,x_o(:,2)*(180/pi))
title('State Response for Observable System')
ylabel('Psi (deg)')
xlabel('Time (s)')
