%%  OBSERVADOR SINCRONIZADO
%% Preamble %%
clc
clear 
close all
rng(12345)
%% Simulation parameters %%
ts = 0.005;
tf = 20;
t = 0:ts:tf;
N = length(t);

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

C = [C1; C2; C3; C4; C5];   % 10×20
OBS = rank(obsv(A, C));

%% Observador %% 
% p=[-145 -150 -155 -160 -165 -170 -175 -180 -185 -150 -155 -160 -165 -170 -175 -180 -170 -180 -190 -195]*0.5;
p = [-2 -3 -4 -5 -6 -7 -8 -9 -10 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13];
L=place(A',C',p)';
KK=L;

%% Condiciones Iniciales %%


% z = [-1;  2; 0; 0;  % robot 1: pos + vel
%      -2; -1; 0; 0;  % robot 2
%       2;  1; 0; 0;  % robot 3
%       2; -2; 0; 0;  % robot 4
%       0;  0; 0; 0]; % robot 5

d = 2.5;  % separación entre líneas (ajusta a tu gusto)

z = [-1.8*d;  0.5; 0; 0;   % robot 1 ROSA
     -1*d;  -1; 0; 0;   % robot 2  verde 
      3;    1; 0; 0;   % robot 3   rojo
      1*d;  -2; 0; 0;   % robot 4  amarillo
      9;  1.5; 0; 0];  % robot 5 morado

% % Initial conditions for observers - diferentes de los reales
% zob1 = z + 0.1*randn(20,1);
% zob2 = z + 0.1*randn(20,1);
% zob3 = z + 0.1*randn(20,1);
% zob4 = z + 0.1*randn(20,1);
% zob5 = z + 0.1*randn(20,1);

zob1 = z + 0.3*randn(20,1);
zob2 = z + 0.3*randn(20,1);
zob3 = z + 0.3*randn(20,1);
zob4 = z + 0.3*randn(20,1);
zob5 = z + 0.3*randn(20,1);

zobtk1 = zob1;
zobtk2 = zob2;
zobtk3 = zob3;
zobtk4 = zob4;
zobtk5 = zob5;

evta1 = 1;
evta2 = 1;
evta3 = 1;
evta4 = 1;
evta5 = 1;

cont1=0;
cont2=0;
cont3=0;
cont4=0;
cont5=0;

ik1=0;
ik2=0;
ik3=0;
ik4=0;
ik5=0;

%% Vectores para guardar %%
Z      = [];
Zob1   = [];
Zob2   = [];
Zob3   = [];
Zob4   = [];
Zob5   = [];

ERR1   = [];
ERR2   = [];
ERR3   = [];
ERR4   = [];
ERR5   = [];
ERRp   = [];

Y1   = [];
Y2   = [];
Y3   = [];
Y4   = [];
Y5   = [];

ZOBTK1 = [];     
ZOBTK2 = [];      
ZOBTK3 = [];      
ZOBTK4 = [];
ZOBTK5 = [];

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

Error1_hist = [];

Zob_p  = [];
Acob_p = [];

%% Tiempos de muestreo
H1 = 0.02; 
H2 = H1*2;  
H3 = H1*3;  
H4 = H1*4;  
H5 = H1*5; 

thetai = 0.1;    % Aumentado para hacer el threshold más restrictivo
bi     = 0.001;     % Aumentado para que eta decaiga más rápido
Gama = 10;

deta1=0;
deta2=0;
deta3=0;
deta4=0;
deta5=0;
%% Ganancias para los términos de enlace %%
Ka = 0.01*eye(20);

%% Observador %%
j=1;
for i=0:ts:tf

    %% EVENT TRIGGERED
    % Cálculo de errores de medición
    dy1=(norm(zob1-zobtk1))^2;
    dy2=(norm(zob2-zobtk2))^2;
    dy3=(norm(zob3-zobtk3))^2;
    dy4=(norm(zob4-zobtk4))^2;
    dy5=(norm(zob5-zobtk5))^2;
    
    % Cálculo de ni_i según la nueva fórmula con r
    ni1 = 0.5 * norm(zobtk1-zobtk1)^2 + ...
          0.5 *  norm(zobtk1-zobtk2)^2 + ...
          0.5 *  norm(zobtk1-zobtk3)^2 + ...
          0.5 * norm(zobtk1-zobtk4)^2 + ...
          0.5 * norm(zobtk1-zobtk5)^2;
    
    ni2 = 0.5 *  norm(zobtk2-zobtk1)^2 + ...
          0.5 * norm(zobtk2-zobtk2)^2 + ...
          0.5 * norm(zobtk2-zobtk3)^2 + ...
          0.5 *  norm(zobtk2-zobtk4)^2 + ...
          0.5 * norm(zobtk2-zobtk5)^2;
    
    ni3 = 0.5 *  norm(zobtk3-zobtk1)^2 + ...
          0.5 *  norm(zobtk3-zobtk2)^2 + ...
          0.5 *  norm(zobtk3-zobtk3)^2 + ...
          0.5 *  norm(zobtk3-zobtk4)^2 + ...
          0.5 *  norm(zobtk3-zobtk5)^2;
    
    ni4 = 0.5 *  norm(zobtk4-zobtk1)^2 + ...
          0.5 *  norm(zobtk4-zobtk2)^2 + ...
          0.5 *  norm(zobtk4-zobtk3)^2 + ...
          0.5 *  norm(zobtk4-zobtk4)^2 + ...
          0.5 *  norm(zobtk4-zobtk5)^2;

    ni5 = 0.5 * norm(zobtk5-zobtk1)^2 + ...
          0.5 * norm(zobtk5-zobtk2)^2 + ...
          0.5 * norm(zobtk5-zobtk3)^2 + ...
          0.5 * norm(zobtk5-zobtk4)^2 + ...
          0.5 * norm(zobtk5-zobtk5)^2;

    y1 = C*z;
    y2 = C*z;
    y3 = C*z;
    y4 = C*z;
    y5 = C*z;


    
    Z  = [Z; z'];
    Zob1 = [Zob1;zob1'];
    Zob2 = [Zob2;zob2'];
    Zob3 = [Zob3;zob3'];
    Zob4 = [Zob4;zob4'];
    Zob5 = [Zob5;zob5'];
    
    vel_deseada = [0.3;0.3];                                        % Velocidades deseadas (cero para posición fija)

    K=place(A(1:4,1:4), B(1:4,1:2), [-2, -3, -4, -5]); % ganancias de control
    u1 = -K * [zeros(2,1); zob1(3:4)  - vel_deseada];
    u2 = -K * [zeros(2,1); zob2(7:8)  - vel_deseada];
    u3 = -K * [zeros(2,1); zob3(11:12)- vel_deseada];
    u4 = -K * [zeros(2,1); zob4(15:16)- vel_deseada];
    u5 = -K * [zeros(2,1); zob5(19:20)- vel_deseada];
    u  = [u1; u2; u3; u4; u5];

    dz = A*z + B*u;


     %% Muestreo 1
        rui1=0.0001;
        ruido1=sqrt(rui1)*randn;
        y1=y1+ruido1;
    if mod(i,H1) == 0
        ikz1    = i;
        ytk1    = y1;
        eta1    = C*zob1 - y1;      % error en el instante de muestreo
    else
        doteta1 = -C*KK*eta1;
        eta1    = eta1 + ts*doteta1;
    end
    Y1=[Y1;y1']; TKZ1=[TKZ1;ikz1]; YTK1=[YTK1;ytk1']; ET1=[ET1;eta1'];

    %% Muestreo 2
        rui2=0.0002;
        ruido2=sqrt(rui2)*randn;
        y2=y2+ruido2;
    if mod(i, H2) == 0
        ikz2    = i;
        ytk2    = y2;
        eta2    = C*zob2 - y2;
    else
        doteta2 = -C*KK*eta2;
        eta2    = eta2 + ts*doteta2;
    end
    Y2=[Y2;y2']; TKZ2=[TKZ2;ikz2]; YTK2=[YTK2;ytk2']; ET2=[ET2;eta2'];

    %% Muestreo 3
        rui3=0.0003;
        ruido3=sqrt(rui3)*randn;
        y3=y3+ruido3;
    if mod(i, H3 ) == 0
        ikz3    = i;
        ytk3    = y3;
        eta3    = C*zob3 - y3;
    else
        doteta3 = -C*KK*eta3;
        eta3    = eta3 + ts*doteta3;
    end
    Y3=[Y3;y3']; TKZ3=[TKZ3;ikz3]; YTK3=[YTK3;ytk3']; ET3=[ET3;eta3'];

    %% Muestreo 4
        rui4=0.0004;
        ruido4=sqrt(rui4)*randn;
        y4=y4+ruido4;
    if mod(i, H4) == 0
        ikz4    = i;
        ytk4    = y4;
        eta4    = C*zob4 - y4;
    else
        doteta4 = -C*KK*eta4;
        eta4    = eta4 + ts*doteta4;
    end
    Y4=[Y4;y4']; TKZ4=[TKZ4;ikz4]; YTK4=[YTK4;ytk4']; ET4=[ET4;eta4'];

    %% Muestreo 5
        rui5=0.0005;
        ruido5=sqrt(rui5)*randn;
        y5=y5+ruido5;
    if mod(i, H5) == 0
        ikz5    = i;
        ytk5    = y5;
        eta5    = C*zob5 - y5;
    else
        doteta5 = -C*KK*eta5;
        eta5    = eta5 + ts*doteta5;
    end
    Y5=[Y5;y5']; TKZ5=[TKZ5;ikz5]; YTK5=[YTK5;ytk5']; ET5=[ET5;eta5'];
        
       %%%%%%%%%%%%%%%%%%%%%%% Comunicación activada por eventos %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Proteger evta_i > 0
    evta1 = max(evta1 + ts*deta1, 1e-10);
    evta2 = max(evta2 + ts*deta2, 1e-10);
    evta3 = max(evta3 + ts*deta3, 1e-10);
    evta4 = max(evta4 + ts*deta4, 1e-10);
    evta5 = max(evta5 + ts*deta5, 1e-10);
    %% Event-Triggered Observer 1
    disp1 = 2*Gama*dy1;
    cond1 = thetai*evta1 + Gama*ni1;

    if disp1 >= cond1
         ik1  = i;
        zobtk1  = zob1;
         ie1  = 1;
       cont1  = cont1+1;
         IE1  = [IE1; ie1'];  
         TK1  = [TK1; ik1];
    end
        ZOBTK1  = [ZOBTK1; zobtk1'];   
    
    %% Event-Triggered Observer 2
    disp2 = 2*Gama*dy2;
    cond2 = thetai*evta2 + Gama*ni2;

    if disp2 >= cond2   
         ik2  = i;
        zobtk2  = zob2;
         ie2  = 2;
       cont2  = cont2+1;
         IE2  = [IE2; ie2'];  
         TK2  = [TK2; ik2];
    end
        ZOBTK2  = [ZOBTK2; zobtk2'];  
    
    %% Event-Triggered Observer 3
    disp3 = 2*Gama*dy3;
    cond3 = thetai*evta3 + Gama*ni3;

    if disp3 >= cond3  
         ik3  = i;
        zobtk3  = zob3;
         ie3  = 3;
       cont3  = cont3+1;
         IE3  = [IE3; ie3'];  
         TK3  = [TK3; ik3];
    end
        ZOBTK3  = [ZOBTK3; zobtk3'];  
    
    %% Event-Triggered Observer 4
    disp4 = 2*Gama*dy4;
    cond4 = thetai*evta4 + Gama*ni4;

    if disp4 >= cond4
         ik4  = i;
        zobtk4  = zob4;
         ie4  = 4;
       cont4  = cont4+1;
         IE4  = [IE4; ie4'];  
         TK4  = [TK4; ik4];
    end
        ZOBTK4  = [ZOBTK4; zobtk4']; 

    %% Event-Triggered Observer 5
    disp5 = 2*Gama*dy5;
    cond5 = thetai*evta5 + Gama*ni5;

    if disp5 >= cond5
         ik5  = i;
        zobtk5  = zob5;
         ie5  = 5;
       cont5  = cont5+1;
         IE5  = [IE5; ie5'];  
         TK5  = [TK5; ik5];
    end
        ZOBTK5  = [ZOBTK5; zobtk5']; 

        %%

        % Dinámica de η_i con comunicación dirigida
    deta1 = -bi*evta1 - 2*Gama*dy1 + Gama*ni1;
    deta2 = -bi*evta2 - 2*Gama*dy2 + Gama*ni2;
    deta3 = -bi*evta3 - 2*Gama*dy3 + Gama*ni3;
    deta4 = -bi*evta4 - 2*Gama*dy4 + Gama*ni4;
    deta5 = -bi*evta5 - 2*Gama*dy5 + Gama*ni5;
     evta1 = evta1 + ts*deta1;
     evta2 = evta2 + ts*deta2;
     evta3 = evta3 + ts*deta3;
     evta4 = evta4 + ts*deta4;
     evta5 = evta5 + ts*deta5;

    % Observadores
    dotZob1 = A*zob1+B*u-KK*eta1+Ka*((zobtk2-zobtk1)+(zobtk3-zobtk1)+(zobtk4-zobtk1)+(zobtk5-zobtk1));
    dotZob2 = A*zob2+B*u-KK*eta2+Ka*((zobtk1-zobtk2)+(zobtk3-zobtk2)+(zobtk4-zobtk2)+(zobtk5-zobtk2));
    dotZob3 = A*zob3+B*u-KK*eta3+Ka*((zobtk1-zobtk3)+(zobtk2-zobtk3)+(zobtk4-zobtk3)+(zobtk5-zobtk3));
    dotZob4 = A*zob4+B*u-KK*eta4+Ka*((zobtk1-zobtk4)+(zobtk2-zobtk4)+(zobtk3-zobtk4)+(zobtk5-zobtk4));
    dotZob5 = A*zob5+B*u-KK*eta5+Ka*((zobtk1-zobtk5)+(zobtk2-zobtk5)+(zobtk3-zobtk5)+(zobtk4-zobtk5));
    
    % --- Integración (Euler) ---
    z    = z    + dz*ts;
    zob1 = zob1 + dotZob1*ts;
    zob2 = zob2 + dotZob2*ts;
    zob3 = zob3 + dotZob3*ts;
    zob4 = zob4 + dotZob4*ts;
    zob5 = zob5 + dotZob5*ts;

    % Error unificado de las estimaciones de posición y velocidad
    
    er1=norm(z - zob1);
    er2=norm(z - zob2);
    er3=norm(z - zob3);
    er4=norm(z - zob4);
    er5=norm(z - zob5);
    
    % Almacenamiento de los errores de posición
    ERR1 = [ERR1; er1'];
    ERR2 = [ERR2; er2'];
    ERR3 = [ERR3; er3'];
    ERR4 = [ERR4; er4'];
    ERR5 = [ERR5; er5'];

    error1 = z - zob1;
    Error1_hist = [Error1_hist; error1'];
    

    %% Estimación promedio de las estimaciones. 

    % Estimación promedio de las estimaciones de posición y velocidad  
        zob_p = (zob1+zob2+zob3+zob4+zob5)/5;
        Zob_p = [Zob_p;zob_p'];

    % Error de la estimación promedio de posición y velocidad
        erp    = norm(z - zob_p);
        ERRp   = [ERRp; erp'];

    % --- Estimación promedio (de los 5) ---
    zob_p = (zob1 + zob2 + zob3 + zob4 + zob5) / 5;
    ERRp(j) = norm(z - zob_p);
 
    j=j+1;
end

%% Extract position data

% Posiciones reales
p1_actual = Z(:,1:2);
p2_actual = Z(:,5:6);
p3_actual = Z(:,9:10);
p4_actual = Z(:,13:14);
p5_actual = Z(:,17:18);

% Estimates from robot 1
p1_est1 = Zob1(:,1:2);
p2_est1 = Zob1(:,5:6);
p3_est1 = Zob1(:,9:10);
p4_est1 = Zob1(:,13:14);
p5_est1 = Zob1(:,17:18);

% Estimates from robot 2
p1_est2 = Zob2(:,1:2);
p2_est2 = Zob2(:,5:6);
p3_est2 = Zob2(:,9:10);
p4_est2 = Zob2(:,13:14);
p5_est2 = Zob2(:,17:18);

% Estimates from robot 3
p1_est3 = Zob3(:,1:2);
p2_est3 = Zob3(:,5:6);
p3_est3 = Zob3(:,9:10);
p4_est3 = Zob3(:,13:14);
p5_est3 = Zob3(:,17:18);

% Estimates from robot 4
p1_est4 = Zob4(:,1:2);
p2_est4 = Zob4(:,5:6);
p3_est4 = Zob4(:,9:10);
p4_est4 = Zob4(:,13:14);
p5_est4 = Zob4(:,17:18);

% Estimates from robot 5
p1_est5 = Zob5(:,1:2);
p2_est5 = Zob5(:,5:6);
p3_est5 = Zob5(:,9:10);
p4_est5 = Zob5(:,13:14);
p5_est5 = Zob5(:,17:18);

%% Create Figure 7.2
figure('Position', [100, 100, 1000, 600])

% Define colors
colors = [0.71 0.43 0.47; % rosa
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


figure;
hold on;
plot(t, ERR1, 'Color', colors(1,:), 'LineWidth', 1.2);
plot(t, ERR2, 'Color', colors(2,:), 'LineWidth', 1.2);
plot(t, ERR3, 'Color', colors(3,:), 'LineWidth', 1.2);
plot(t, ERR4, 'Color', colors(4,:), 'LineWidth', 1.2);
plot(t, ERR5, 'Color', colors(5,:), 'LineWidth', 1.2);
hold off;
xlabel('Tiempo (s)');
ylabel('Error de estimación (norma)');
grid on;
title('Evolución del error de estimación (medición global)');
legend('Robot 1','Robot 2','Robot 3','Robot 4','Robot 5');
%%
figure('Units', 'centimeters', 'Position', [0 0 17 12]);
    plot(t, Z(:,1), 'k', 'LineWidth', 2); hold on
    plot(t, Zob1(:,1), '--', 'Color', colors(1,:), 'LineWidth', 2);
    plot(t, Zob2(:,1), '--', 'Color', colors(2,:), 'LineWidth', 2);
    plot(t, Zob3(:,1), '--', 'Color', colors(3,:), 'LineWidth', 2);
    plot(t, Zob4(:,1), '--', 'Color', colors(4,:), 'LineWidth', 2);
    plot(t, Zob5(:,1), '--', 'Color', colors(5,:), 'LineWidth', 2);
    grid on
    ylabel('Y', 'FontSize', 11, 'Interpreter','latex');
    xlabel('Time [s]','FontSize', 11, 'Interpreter','latex');
    legend('$x$','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$','$\hat{x}_5$','FontSize', 11, 'Interpreter','latex'); % Ajusta los nombres según corresponda
    set(gca, 'FontSize', 11, 'FontName', 'Times', 'LineWidth', 0.5)

    %% 
% ╔═══════════════════════════════════════════════════════════════════════╗
% ║                         CÁLCULO DE ERRORES                            ║
% ╚═══════════════════════════════════════════════════════════════════════╝
error_obs1 = Zob1(:,:) - Z(:,:);  % Error observador 1
error_obs2 = Zob2(:,:) - Z(:,:);  % Error observador 2
error_obs3 = Zob3(:,:) - Z(:,:);  % Error observador 3
error_obs4 = Zob4(:,:) - Z(:,:);  % Error observador 4
error_obs5 = Zob5(:,:) - Z(:,:);  % Error observador 5


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