e%% Preamble %%
clc
clear 
close all

%% Simulation parameters %%
ts = 0.001;
tf = 30;
t = 0:ts:tf;

%% System matrices A, B
A = kron(eye(5), [zeros(2,2), eye(2); zeros(2,4)]);

B_block = [zeros(2,2); eye(2)];
B = kron(eye(5), B_block);

S1 = [1,  0,  0,  0,  0];
S2 = [1, -1,  0,  0,  0];
S3 = [1,  0, -1,  0,  0];
S4 = [0,  1,  0, -1,  0];
S5 = [0,  0,  1,  0, -1];

% Medir solo POSICIÓN (primeros 2 estados de cada robot)
C1 = [kron(S1, [1 0 0 0; 0 1 0 0])];  % 2x20
C2 = [kron(S2, [1 0 0 0; 0 1 0 0])];
C3 = [kron(S3, [1 0 0 0; 0 1 0 0])];
C4 = [kron(S4, [1 0 0 0; 0 1 0 0])];
C5 = [kron(S5, [1 0 0 0; 0 1 0 0])];
%% Communication topology %%

La = [0 1 1 0 0
      1 0 0 1 0
      1 0 0 0 1
      0 1 0 0 0
      0 0 1 0 0];

%% Transformation matrices for observability decomposition
[~,~,~,T1,k1] = obsvf(rot90(A,2),0,rot90(C1,2)); 
[~,~,~,T2,k2] = obsvf(rot90(A,2),0,rot90(C2,2)); 
[~,~,~,T3,k3] = obsvf(rot90(A,2),0,rot90(C3,2)); 
[~,~,~,T4,k4] = obsvf(rot90(A,2),0,rot90(C4,2)); 
[~,~,~,T5,k5] = obsvf(rot90(A,2),0,rot90(C5,2)); 

T1 = rot90(T1,2)';
T2 = rot90(T2,2)';
T3 = rot90(T3,2)';
T4 = rot90(T4,2)';
T5 = rot90(T5,2)';

v1o = rank(obsv(A,C1)); % tamaño 4
v2o = rank(obsv(A,C2)); % tamaño 4
v3o = rank(obsv(A,C3)); % tamaño 4
v4o = rank(obsv(A,C4)); % tamaño 4
v5o = rank(obsv(A,C5)); % tamaño 4

n = size(A,1); 

Abar1 = T1'*A*T1; 
Abar2 = T2'*A*T2; 
Abar3 = T3'*A*T3; 
Abar4 = T4'*A*T4;
Abar5 = T5'*A*T5;

A1o = Abar1(1:v1o,1:v1o); 
A2o = Abar2(1:v2o,1:v2o); 
A3o = Abar3(1:v3o,1:v3o); 
A4o = Abar4(1:v4o,1:v4o); 
A5o = Abar5(1:v5o,1:v5o); 

Cbar1 = C1*T1; 
Cbar2 = C2*T2; 
Cbar3 = C3*T3; 
Cbar4 = C4*T4; 
Cbar5 = C5*T5; 

C1o = Cbar1(:,1:v1o);
C2o = Cbar2(:,1:v2o);
C3o = Cbar3(:,1:v3o);
C4o = Cbar4(:,1:v4o);
C5o = Cbar5(:,1:v5o);

%% Pole placement - todos necesitan 4 polos
polos_L1o = [-15,   -14,   -13,   -12];
polos_L2o = [-15,   -14,   -13,   -12];
polos_L3o = [-16,   -15,   -14,   -13];   
polos_L4o = [-16,   -15,   -14,   -13];
polos_L5o = [-16,   -15,   -14,   -13];   

% Calculate gains L for observers and M for communication
L1o = place(A1o', C1o', polos_L1o)';
L2o = place(A2o', C2o', polos_L2o)';
L3o = place(A3o', C3o', polos_L3o)';
L4o = place(A4o', C4o', polos_L4o)';
L5o = place(A5o', C5o', polos_L5o)';
 
L1 = T1 * [L1o; zeros(16, 2)];
L2 = T2 * [L2o; zeros(16, 2)];
L3 = T3 * [L3o; zeros(16, 2)];
L4 = T4 * [L4o; zeros(16, 2)];
L5 = T5 * [L5o; zeros(16, 2)];

M1 = T1 * [zeros(v1o,v1o) zeros(v1o,n-v1o); zeros(n-v1o,v1o) eye(n-v1o)] * T1';
M2 = T2 * [zeros(v2o,v2o) zeros(v2o,n-v2o); zeros(n-v2o,v2o) eye(n-v2o)] * T2';
M3 = T3 * [zeros(v3o,v3o) zeros(v3o,n-v3o); zeros(n-v3o,v3o) eye(n-v3o)] * T3';
M4 = T4 * [zeros(v4o,v4o) zeros(v4o,n-v4o); zeros(n-v4o,v4o) eye(n-v4o)] * T4';
M5 = T5 * [zeros(v5o,v5o) zeros(v5o,n-v5o); zeros(n-v5o,v5o) eye(n-v5o)] * T5';

x = [-1;  2; 0; 0;  % robot 1: pos + vel
     -2; -1; 0; 0;  % robot 2
      2;  1; 0; 0;  % robot 3
      2; -2; 0; 0;  % robot 4
      0;  0; 0; 0]; % robot 5

%% Initial conditions for observers - diferentes de los reales
x1 = x + 0.1*randn(20,1);
x2 = x + 0.1*randn(20,1);
x3 = x + 0.1*randn(20,1);
x4 = x + 0.1*randn(20,1);
x5 = x + 0.1*randn(20,1);
%% Storage vectors
X   = [];
X1  = [];   
X2  = [];   
X3  = [];   
X4  = [];
X5  = [];

Y1   = [];
Y2   = [];
Y3   = [];
Y4   = [];
Y5   = [];

XTK1 = [];     
XTK2 = [];      
XTK3 = [];      
XTK4 = [];
XTK5 = [];

TK1  = []; 
TK2  = []; 
TK3  = []; 
TK4  = []; 
TK5  = [];

IE1  = []; 
IE2  = [];
IE3  = []; 
IE4  = []; 
IE5  = [];

YTK1 = [];
YTK2 = [];
YTK3 = [];
YTK4 = [];
YTK5 = [];

TKZ1 = [];
TKZ2 = [];
TKZ3 = [];
TKZ4 = [];
TKZ5 = [];

ET1 = [];
ET2 = [];
ET3 = [];
ET4 = [];
ET5 = [];

Gama = 35; % Ganancia de acoplamiento

%% Tiempos de muestreo → convertir a pasos enteros
H1 = 0.04; 
H2 = H1*2;  
H3 = H1*3;  
H4 = H1*4;  
H5 = H1*5; 
%% ── Inicialización de eta ANTES del loop ──────────────────────────────
eta1 = C1*x1 - C1*x;
eta2 = C2*x2 - C2*x;
eta3 = C3*x3 - C3*x;
eta4 = C4*x4 - C4*x;
eta5 = C5*x5 - C5*x;

ikz1=0; ikz2=0; ikz3=0; ikz4=0; ikz5=0;
ytk1=C1*x; ytk2=C2*x; ytk3=C3*x; ytk4=C4*x; ytk5=C5*x;

%% Simulation loop

for i = 0:ts:tf
    X  = [X; x'];
    X1 = [X1; x1'];
    X2 = [X2; x2'];
    X3 = [X3; x3'];
    X4 = [X4; x4'];
    X5 = [X5; x5'];
    
     %% Controlador
    vel_deseada = [0.3;0.3];                                        % Velocidades deseadas (cero para posición fija)

    K=place(A(1:4,1:4), B(1:4,1:2), [-2, -3, -4, -5]); % ganancias de control
    u1 = -K * [zeros(2,1); x1(3:4)  - vel_deseada];
    u2 = -K * [zeros(2,1); x2(7:8)  - vel_deseada];
    u3 = -K * [zeros(2,1); x3(11:12)- vel_deseada];
    u4 = -K * [zeros(2,1); x4(15:16)- vel_deseada];
    u5 = -K * [zeros(2,1); x5(19:20)- vel_deseada];
    u  = [u1; u2; u3; u4; u5];

    dx = A*x + B*u;

    y1 = C1*x;
    y2 = C2*x;
    y3 = C3*x;
    y4 = C4*x;
    y5 = C5*x;

    y1hat = C1*x1;
    y2hat = C2*x2;
    y3hat = C3*x3;
    y4hat = C4*x4;
    y5hat = C5*x5;

    %% Muestreo 1
    if mod(i,H1) == 0
        ikz1    = i;
        ytk1    = y1;
        eta1    = C1*x1 - y1;      % error en el instante de muestreo
    else
        doteta1 = -C1*L1*eta1;
        eta1    = eta1 + ts*doteta1;
    end
    Y1=[Y1;y1']; TKZ1=[TKZ1;ikz1]; YTK1=[YTK1;ytk1']; ET1=[ET1;eta1'];

    %% Muestreo 2
    if mod(i, H2) == 0
        ikz2    = i;
        ytk2    = y2;
        eta2    = C2*x2 - y2;
    else
        doteta2 = -C2*L2*eta2;
        eta2    = eta2 + ts*doteta2;
    end
    Y2=[Y2;y2']; TKZ2=[TKZ2;ikz2]; YTK2=[YTK2;ytk2']; ET2=[ET2;eta2'];

    %% Muestreo 3
    if mod(i, H3 ) == 0
        ikz3    = i;
        ytk3    = y3;
        eta3    = C3*x3 - y3;
    else
        doteta3 = -C3*L3*eta3;
        eta3    = eta3 + ts*doteta3;
    end
    Y3=[Y3;y3']; TKZ3=[TKZ3;ikz3]; YTK3=[YTK3;ytk3']; ET3=[ET3;eta3'];

    %% Muestreo 4
    if mod(i, H4) == 0
        ikz4    = i;
        ytk4    = y4;
        eta4    = C4*x4 - y4;
    else
        doteta4 = -C4*L4*eta4;
        eta4    = eta4 + ts*doteta4;
    end
    Y4=[Y4;y4']; TKZ4=[TKZ4;ikz4]; YTK4=[YTK4;ytk4']; ET4=[ET4;eta4'];

    %% Muestreo 5
    if mod(i, H5) == 0
        ikz5    = i;
        ytk5    = y5;
        eta5    = C5*x5 - y5;
    else
        doteta5 = -C5*L5*eta5;
        eta5    = eta5 + ts*doteta5;
    end
    Y5=[Y5;y5']; TKZ5=[TKZ5;ikz5]; YTK5=[YTK5;ytk5']; ET5=[ET5;eta5'];


    dx1 = A*x1 + B*u - L1*eta1 + ...
          Gama*M1*(La(1,1)*(x1-x1) + La(1,2)*(x2-x1) + La(1,3)*(x3-x1) + La(1,4)*(x4-x1) + La(1,5)*(x5-x1));
    dx2 = A*x2 + B*u - L2*eta2 + ...
          Gama*M2*(La(2,1)*(x1-x2) + La(2,2)*(x2-x2) + La(2,3)*(x3-x2) + La(2,4)*(x4-x2) + La(2,5)*(x5-x2));
    dx3 = A*x3 + B*u - L3*eta3 + ...
          Gama*M3*(La(3,1)*(x1-x3) + La(3,2)*(x2-x3) + La(3,3)*(x3-x3) + La(3,4)*(x4-x3) + La(3,5)*(x5-x3));
    dx4 = A*x4 + B*u - L4*eta4 + ...
          Gama*M4*(La(4,1)*(x1-x4) + La(4,2)*(x2-x4) + La(4,3)*(x3-x4) + La(4,4)*(x4-x4) + La(4,5)*(x5-x4));
    dx5 = A*x5 + B*u - L5*eta5 + ...
          Gama*M5*(La(5,1)*(x1-x5) + La(5,2)*(x2-x5) + La(5,3)*(x3-x5) + La(5,4)*(x4-x5) + La(5,5)*(x5-x5));
    
    % Numerical integration (Euler)
    x  = x  + ts*dx;
    x1 = x1 + ts*dx1;
    x2 = x2 + ts*dx2;
    x3 = x3 + ts*dx3;
    x4 = x4 + ts*dx4;
    x5 = x5 + ts*dx5;
end  

%% Extract position data
% Actual positions
p1_actual = X(:,1:2);
p2_actual = X(:,5:6);
p3_actual = X(:,9:10);
p4_actual = X(:,13:14);
p5_actual = X(:,17:18);

% Estimates from robot 1
p1_est1 = X1(:,1:2);
p2_est1 = X1(:,5:6);
p3_est1 = X1(:,9:10);
p4_est1 = X1(:,13:14);
p5_est1 = X1(:,17:18);

% Estimates from robot 2
p1_est2 = X2(:,1:2);
p2_est2 = X2(:,5:6);
p3_est2 = X2(:,9:10);
p4_est2 = X2(:,13:14);
p5_est2 = X2(:,17:18);

% Estimates from robot 3
p1_est3 = X3(:,1:2);
p2_est3 = X3(:,5:6);
p3_est3 = X3(:,9:10);
p4_est3 = X3(:,13:14);
p5_est3 = X3(:,17:18);

% Estimates from robot 4
p1_est4 = X4(:,1:2);
p2_est4 = X4(:,5:6);
p3_est4 = X4(:,9:10);
p4_est4 = X4(:,13:14);
p5_est4 = X4(:,17:18);

% Estimates from robot 5
p1_est5 = X5(:,1:2);
p2_est5 = X5(:,5:6);
p3_est5 = X5(:,9:10);
p4_est5 = X5(:,13:14);
p5_est5 = X5(:,17:18);

%% Create Figure 7.2
figure('Position', [100, 100, 1000, 600])

% Define colors
colors = [0.71 0.43 0.47;       % rosa
          0.46 0.67 0.19; % verde
          0.85 0.33 0.10; % morado
          0.93 0.69 0.13; % naranja
          0.49 0.18 0.56];% amarillo

% Subplot (a): Estimate from robot 1
subplot(2, 3, 1)
hold on
% Actual trajectories (black dashed)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p2_actual(:,1), p2_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p3_actual(:,1), p3_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p4_actual(:,1), p4_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p5_actual(:,1), p5_actual(:,2), 'k--', 'LineWidth', 1.5)

% Estimated trajectories from robot 1 (colored solid)
plot(p1_est1(:,1), p1_est1(:,2), 'Color', colors(1,:), 'LineWidth', 2)
plot(p2_est1(:,1), p2_est1(:,2), 'Color', colors(2,:), 'LineWidth', 2)
plot(p3_est1(:,1), p3_est1(:,2), 'Color', colors(3,:), 'LineWidth', 2)
plot(p4_est1(:,1), p4_est1(:,2), 'Color', colors(4,:), 'LineWidth', 2)
plot(p5_est1(:,1), p5_est1(:,2), 'Color', colors(5,:), 'LineWidth', 2)

hold off
xlabel('x')
ylabel('y')
title('(a) Estimate from robot 1')
axis equal
grid on
xlim([-5 15])
ylim([-5 15])

% Subplot (b): Estimate from robot 2
subplot(2, 3, 2)
hold on
% Actual trajectories (black dashed)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p2_actual(:,1), p2_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p3_actual(:,1), p3_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p4_actual(:,1), p4_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p5_actual(:,1), p5_actual(:,2), 'k--', 'LineWidth', 1.5)

% Estimated trajectories from robot 2 (colored solid)
plot(p1_est2(:,1), p1_est2(:,2), 'Color', colors(1,:), 'LineWidth', 2)
plot(p2_est2(:,1), p2_est2(:,2), 'Color', colors(2,:), 'LineWidth', 2)
plot(p3_est2(:,1), p3_est2(:,2), 'Color', colors(3,:), 'LineWidth', 2)
plot(p4_est2(:,1), p4_est2(:,2), 'Color', colors(4,:), 'LineWidth', 2)
plot(p5_est2(:,1), p5_est2(:,2), 'Color', colors(5,:), 'LineWidth', 2)

hold off
xlabel('x')
ylabel('y')
title('(b) Estimate from robot 2')
axis equal
grid on
xlim([-5 15])
ylim([-5 15])

% Subplot (c): Estimate from robot 3
subplot(2, 3, 3)
hold on
% Actual trajectories (black dashed)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p2_actual(:,1), p2_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p3_actual(:,1), p3_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p4_actual(:,1), p4_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p5_actual(:,1), p5_actual(:,2), 'k--', 'LineWidth', 1.5)

% Estimated trajectories from robot 3 (colored solid)
plot(p1_est3(:,1), p1_est3(:,2), 'Color', colors(1,:), 'LineWidth', 2)
plot(p2_est3(:,1), p2_est3(:,2), 'Color', colors(2,:), 'LineWidth', 2)
plot(p3_est3(:,1), p3_est3(:,2), 'Color', colors(3,:), 'LineWidth', 2)
plot(p4_est3(:,1), p4_est3(:,2), 'Color', colors(4,:), 'LineWidth', 2)
plot(p5_est3(:,1), p5_est3(:,2), 'Color', colors(5,:), 'LineWidth', 2)

hold off
xlabel('x')
ylabel('y')
title('(c) Estimate from robot 3')
axis equal
grid on
xlim([-5 15])
ylim([-5 15])

% Subplot (d): Estimate from robot 4
subplot(2, 3, 4)
hold on
% Actual trajectories (black dashed)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p2_actual(:,1), p2_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p3_actual(:,1), p3_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p4_actual(:,1), p4_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p5_actual(:,1), p5_actual(:,2), 'k--', 'LineWidth', 1.5)

% Estimated trajectories from robot 4 (colored solid)
plot(p1_est4(:,1), p1_est4(:,2), 'Color', colors(1,:), 'LineWidth', 2)
plot(p2_est4(:,1), p2_est4(:,2), 'Color', colors(2,:), 'LineWidth', 2)
plot(p3_est4(:,1), p3_est4(:,2), 'Color', colors(3,:), 'LineWidth', 2)
plot(p4_est4(:,1), p4_est4(:,2), 'Color', colors(4,:), 'LineWidth', 2)
plot(p5_est4(:,1), p5_est4(:,2), 'Color', colors(5,:), 'LineWidth', 2)

hold off
xlabel('x')
ylabel('y')
title('(d) Estimate from robot 4')
axis equal
grid on
xlim([-5 15])
ylim([-5 15])

% Subplot (e): Estimate from robot 5
subplot(2, 3, 5)
hold on
% Actual trajectories (black dashed)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p2_actual(:,1), p2_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p3_actual(:,1), p3_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p4_actual(:,1), p4_actual(:,2), 'k--', 'LineWidth', 1.5)
plot(p5_actual(:,1), p5_actual(:,2), 'k--', 'LineWidth', 1.5)

% Estimated trajectories from robot 5 (colored solid)
plot(p1_est5(:,1), p1_est5(:,2), 'Color', colors(1,:), 'LineWidth', 2)
plot(p2_est5(:,1), p2_est5(:,2), 'Color', colors(2,:), 'LineWidth', 2)
plot(p3_est5(:,1), p3_est5(:,2), 'Color', colors(3,:), 'LineWidth', 2)
plot(p4_est5(:,1), p4_est5(:,2), 'Color', colors(4,:), 'LineWidth', 2)
plot(p5_est5(:,1), p5_est5(:,2), 'Color', colors(5,:), 'LineWidth', 2)

hold off
xlabel('x')
ylabel('y')
title('(e) Estimate from robot 5')
axis equal
grid on
xlim([-5 15])
ylim([-5 15])
%% 
% ╔═══════════════════════════════════════════════════════════════════════╗
% ║                         CÁLCULO DE ERRORES                            ║
% ╚═══════════════════════════════════════════════════════════════════════╝
error_obs1 = X1(:,:) - X(:,:);  % Error observador 1
error_obs2 = X2(:,:) - X(:,:);  % Error observador 2
error_obs3 = X3(:,:) - X(:,:);  % Error observador 3
error_obs4 = X4(:,:) - X(:,:);  % Error observador 4
error_obs5 = X5(:,:) - X(:,:);  % Error observador 5


WE11= vecnorm([error_obs1(:,1:2)],  2, 2); %observador 1 del robot 1
WE21= vecnorm([error_obs2(:,1:2)],  2, 2); %observador 2 del robot 1
WE31= vecnorm([error_obs3(:,1:2)],  2, 2); %observador 3 del robot 1
WE41= vecnorm([error_obs4(:,1:2)],  2, 2); %observador 4 del robot 1
WE51= vecnorm([error_obs5(:,1:2)],  2, 2); %observador 5 del robot 1

WE12= vecnorm([error_obs1(:,5:6)],  2, 2); %observador 1 del robot 2
WE22= vecnorm([error_obs2(:,5:6)],  2, 2); %observador 2 del robot 2
WE32= vecnorm([error_obs3(:,5:6)],  2, 2); %observador 3 del robot 2
WE42= vecnorm([error_obs4(:,5:6)],  2, 2); %observador 4 del robot 2
WE52= vecnorm([error_obs5(:,5:6)],  2, 2); %observador 5 del robot 2

WE13= vecnorm([error_obs1(:,9:10)],  2, 2); %observador 1 del robot 3
WE23= vecnorm([error_obs2(:,9:10)],  2, 2); %observador 2 del robot 3
WE33= vecnorm([error_obs3(:,9:10)],  2, 2); %observador 3 del robot 3
WE43= vecnorm([error_obs4(:,9:10)],  2, 2); %observador 4 del robot 3
WE53= vecnorm([error_obs5(:,9:10)],  2, 2); %observador 5 del robot 3

WE14= vecnorm([error_obs1(:,13:14)],  2, 2); %observador 1 del robot 4
WE24= vecnorm([error_obs2(:,13:14)],  2, 2); %observador 2 del robot 4
WE34= vecnorm([error_obs3(:,13:14)],  2, 2); %observador 3 del robot 4
WE44= vecnorm([error_obs4(:,13:14)],  2, 2); %observador 4 del robot 4
WE54= vecnorm([error_obs5(:,13:14)],  2, 2); %observador 5 del robot 4

WE15= vecnorm([error_obs1(:,17:18)],  2, 2); %observador 1 del robot 5
WE25= vecnorm([error_obs2(:,17:18)],  2, 2); %observador 2 del robot 5
WE35= vecnorm([error_obs3(:,17:18)],  2, 2); %observador 3 del robot 5
WE45= vecnorm([error_obs4(:,17:18)],  2, 2); %observador 4 del robot 5
WE55= vecnorm([error_obs5(:,17:18)],  2, 2); %observador 5 del robot 5

dt       = ts;
N_robots = 5;
N_obs    = 5;

% Inicializar matrices de métricas
IAE_pos  = zeros(N_obs, N_robots);
ITAE_pos = zeros(N_obs, N_robots);
ISE_pos  = zeros(N_obs, N_robots);

% Llenar para cada observador y robot
% Observador 1
IAE_pos(1,1)  = sum(WE11) * dt;
ITAE_pos(1,1) = sum(t .* WE11') * dt;
ISE_pos(1,1)  = sum(WE11.^2) * dt;

IAE_pos(1,2)  = sum(WE12) * dt;
ITAE_pos(1,2) = sum(t .* WE12') * dt;
ISE_pos(1,2)  = sum(WE12.^2) * dt;

IAE_pos(1,3)  = sum(WE13) * dt;
ITAE_pos(1,3) = sum(t .* WE13') * dt;
ISE_pos(1,3)  = sum(WE13.^2) * dt;

IAE_pos(1,4)  = sum(WE14) * dt;
ITAE_pos(1,4) = sum(t .* WE14') * dt;
ISE_pos(1,4)  = sum(WE14.^2) * dt;

IAE_pos(1,5)  = sum(WE15) * dt;
ITAE_pos(1,5) = sum(t .* WE15') * dt;
ISE_pos(1,5)  = sum(WE15.^2) * dt;

% Observador 2
IAE_pos(2,1)  = sum(WE21) * dt;
ITAE_pos(2,1) = sum(t .* WE21') * dt;
ISE_pos(2,1)  = sum(WE21.^2) * dt;

IAE_pos(2,2)  = sum(WE22) * dt;
ITAE_pos(2,2) = sum(t .* WE22') * dt;
ISE_pos(2,2)  = sum(WE22.^2) * dt;

IAE_pos(2,3)  = sum(WE23) * dt;
ITAE_pos(2,3) = sum(t .* WE23') * dt;
ISE_pos(2,3)  = sum(WE23.^2) * dt;

IAE_pos(2,4)  = sum(WE24) * dt;
ITAE_pos(2,4) = sum(t .* WE24') * dt;
ISE_pos(2,4)  = sum(WE24.^2) * dt;

IAE_pos(2,5)  = sum(WE25) * dt;
ITAE_pos(2,5) = sum(t .* WE25') * dt;
ISE_pos(2,5)  = sum(WE25.^2) * dt;

% Observador 3
IAE_pos(3,1)  = sum(WE31) * dt;
ITAE_pos(3,1) = sum(t .* WE31') * dt;
ISE_pos(3,1)  = sum(WE31.^2) * dt;

IAE_pos(3,2)  = sum(WE32) * dt;
ITAE_pos(3,2) = sum(t .* WE32') * dt;
ISE_pos(3,2)  = sum(WE32.^2) * dt;

IAE_pos(3,3)  = sum(WE33) * dt;
ITAE_pos(3,3) = sum(t .* WE33') * dt;
ISE_pos(3,3)  = sum(WE33.^2) * dt;

IAE_pos(3,4)  = sum(WE34) * dt;
ITAE_pos(3,4) = sum(t .* WE34') * dt;
ISE_pos(3,4)  = sum(WE34.^2) * dt;

IAE_pos(3,5)  = sum(WE35) * dt;
ITAE_pos(3,5) = sum(t .* WE35') * dt;
ISE_pos(3,5)  = sum(WE35.^2) * dt;

% Observador 4
IAE_pos(4,1)  = sum(WE41) * dt;
ITAE_pos(4,1) = sum(t .* WE41') * dt;
ISE_pos(4,1)  = sum(WE41.^2) * dt;

IAE_pos(4,2)  = sum(WE42) * dt;
ITAE_pos(4,2) = sum(t .* WE42') * dt;
ISE_pos(4,2)  = sum(WE42.^2) * dt;

IAE_pos(4,3)  = sum(WE43) * dt;
ITAE_pos(4,3) = sum(t .* WE43') * dt;
ISE_pos(4,3)  = sum(WE43.^2) * dt;

IAE_pos(4,4)  = sum(WE44) * dt;
ITAE_pos(4,4) = sum(t .* WE44') * dt;
ISE_pos(4,4)  = sum(WE44.^2) * dt;

IAE_pos(4,5)  = sum(WE45) * dt;
ITAE_pos(4,5) = sum(t .* WE45') * dt;
ISE_pos(4,5)  = sum(WE45.^2) * dt;

% Observador 5
IAE_pos(5,1)  = sum(WE51) * dt;
ITAE_pos(5,1) = sum(t .* WE51') * dt;
ISE_pos(5,1)  = sum(WE51.^2) * dt;

IAE_pos(5,2)  = sum(WE52) * dt;
ITAE_pos(5,2) = sum(t .* WE52') * dt;
ISE_pos(5,2)  = sum(WE52.^2) * dt;

IAE_pos(5,3)  = sum(WE53) * dt;
ITAE_pos(5,3) = sum(t .* WE53') * dt;
ISE_pos(5,3)  = sum(WE53.^2) * dt;

IAE_pos(5,4)  = sum(WE54) * dt;
ITAE_pos(5,4) = sum(t .* WE54') * dt;
ISE_pos(5,4)  = sum(WE54.^2) * dt;

IAE_pos(5,5)  = sum(WE55) * dt;
ITAE_pos(5,5) = sum(t .* WE55') * dt;
ISE_pos(5,5)  = sum(WE55.^2) * dt;

% Luego puedes promediar por observador (sobre los robotes) o viceversa
promedio_IAE_obs = mean(IAE_pos, 2); % promedio por observador
promedio_IAE_robot = mean(IAE_pos, 1); % promedio por robot

% Luego puedes promediar por observador (sobre los robotes) o viceversa
promedio_ITAE_obs = mean(ITAE_pos, 2); % promedio por observador
promedio_ITAE_robot = mean(ITAE_pos, 1); % promedio por robot

% Luego puedes promediar por observador (sobre los robotes) o viceversa
promedio_ISE_obs = mean(ISE_pos, 2); % promedio por observador
promedio_ISE_robot = mean(ISE_pos, 1); % promedio por robot

%%  TABLA CON ERRORES
tabla_metricas = table(...
    round(promedio_ITAE_obs, 5), ...
    round(promedio_IAE_obs, 5), ...
    round(promedio_ISE_obs, 5), ...
    'VariableNames', {'ITAE Average', 'IAE Average', 'ISE Average'}, ...
    'RowNames', {'x1_hat', 'x2_hat', 'x3_hat', 'x4_hat','x5_hat'}); 
disp('=== PROMEDIO DE MÉTRICAS ===')
disp(tabla_metricas)
%%
figure('Units', 'centimeters', 'Position', [0 0 17 12]);
    plot(t, X(:,1), 'k', 'LineWidth', 2); hold on
    plot(t, X1(:,1), '--', 'Color', colors(1,:), 'LineWidth', 2);
    plot(t, X2(:,1), '--', 'Color', colors(2,:), 'LineWidth', 2);
    plot(t, X3(:,1), '--', 'Color', colors(3,:), 'LineWidth', 2);
    plot(t, X4(:,1), '--', 'Color', colors(4,:), 'LineWidth', 2);
    plot(t, X5(:,1), '--', 'Color', colors(5,:), 'LineWidth', 2);
    grid on
    ylabel('Y', 'FontSize', 11, 'Interpreter','latex');
    xlabel('Time [s]','FontSize', 11, 'Interpreter','latex');
    legend('$x$','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$','$\hat{x}_5$','FontSize', 11, 'Interpreter','latex'); % Ajusta los nombres según corresponda
    set(gca, 'FontSize', 11, 'FontName', 'Times', 'LineWidth', 0.5)
