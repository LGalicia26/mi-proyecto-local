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

%% Matrices para un solo sistema (1 robot/estructura)
% El bloque base representa un doble integrador en 2D (4 estados en total)
% Estados: [pos_x; pos_y; vel_x; vel_y]

A = [zeros(2,2),  eye(2); 
     zeros(2,2),  zeros(2,2)]; % Matriz 4x4

B = [zeros(2,2); 
       eye(2)];                  % Matriz 4x2

% Medir solo POSICIÓN (los primeros 2 estados)
C = [eye(2), zeros(2,2)];      % Matriz 2x4

%% Verificación de Observabilidad
OBS_single = rank(obsv(A, C)); % 4

%% Observador %% 
p = [-2 -3 -4 -5];
KK=place(A',C',p)';

z = [2.5;  0.5; 0; 0];   % robot 1 ROSA
 
% Initial conditions for observers - diferentes de los reales
zob1 = z + 0.3*randn(4,1);
zob2 = z + 0.3*randn(4,1);
zob3 = z + 0.3*randn(4,1);


zobtk1 = zob1;
zobtk2 = zob2;
zobtk3 = zob3;

evta1 = 1;
evta2 = 1;
evta3 = 1;


cont1=0;
cont2=0;
cont3=0;


ik1=0;
ik2=0;
ik3=0;


%% Vectores para guardar %%
Z      = [];
Zob1   = [];
Zob2   = [];
Zob3   = [];


ERR1   = [];
ERR2   = [];
ERR3   = [];

ERRp   = [];

Y1   = [];
Y2   = [];
Y3   = [];


ZOBTK1 = [];     
ZOBTK2 = [];      
ZOBTK3 = [];      


TK1  = []; 
TK2  = []; 
TK3  = []; 

IE1  = []; 
IE2  = [];
IE3  = []; 

YTK1 = [];
YTK2 = [];
YTK3 = [];

TKZ1 = [];
TKZ2 = [];
TKZ3 = [];

ET1 = [];
ET2 = [];
ET3 = [];


Error1_hist = [];

Zob_p  = [];
Acob_p = [];

%% Tiempos de muestreo
H1 = 0.02; 
H2 = H1*2;  
H3 = H1*3;  
 
deta1=0;
deta2=0;
deta3=0;

%% Ganancias para los términos de enlace %%
Ka = 0.01*eye(4);

%% Observador %%
j=1;
for i=0:ts:tf

      
    y1 = C*z;
    y2 = C*z;
    y3 = C*z;
    
    Z  = [Z; z'];
    Zob1 = [Zob1;zob1'];
    Zob2 = [Zob2;zob2'];
    Zob3 = [Zob3;zob3'];
    

    vel_deseada = [0.3;0.3];                                        % Velocidades deseadas (cero para posición fija)

    K=place(A(1:4,1:4), B(1:4,1:2), [-2, -3, -4, -5]); % ganancias de control
    u1 = -K * [zeros(2,1); zob1(3:4)  - vel_deseada];
      
    u  = u1;

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

   


    % Observadores
    dotZob1 = A*zob1+B*u-KK*eta1+Ka*((zob2-zob1)+(zob3-zob1));
    dotZob2 = A*zob2+B*u-KK*eta2+Ka*((zob1-zob2)+(zob3-zob2));
    dotZob3 = A*zob3+B*u-KK*eta3+Ka*((zob1-zob3)+(zob2-zob3));
    
    % --- Integración (Euler) ---
    z    = z    + dz*ts;
    zob1 = zob1 + dotZob1*ts;
    zob2 = zob2 + dotZob2*ts;
    zob3 = zob3 + dotZob3*ts;
    

    % Error unificado de las estimaciones de posición y velocidad

    er1=norm(z - zob1);
    er2=norm(z - zob2);
    er3=norm(z - zob3);
    

    % Almacenamiento de los errores de posición
    ERR1 = [ERR1; er1'];
    ERR2 = [ERR2; er2'];
    ERR3 = [ERR3; er3'];
    
    error1 = z - zob1;
    Error1_hist = [Error1_hist; error1'];


    %% Estimación promedio de las estimaciones. 

    % Estimación promedio de las estimaciones de posición y velocidad  
        zob_p = (zob1+zob2+zob3)/3;
        Zob_p = [Zob_p;zob_p'];

    % Error de la estimación promedio de posición y velocidad
        erp    = norm(z - zob_p);
        ERRp   = [ERRp; erp'];

        p1_avg = Zob_p(:, 1:2);

    j=j+1;
end

%% Extract position data

% Posiciones reales
p1_actual = Z(:,1:2);


% Estimates from robot 1
p1_est1 = Zob1(:,1:2);


% Estimates from robot 2
p1_est2 = Zob2(:,1:2);


% Estimates from robot 3
p1_est3 = Zob3(:,1:2);


%% Create Figure 7.2
figure('Position', [100, 100, 1000, 600])

% Define colors
colors = [0.71 0.43 0.47; % rosa
          0.46 0.67 0.19; % verde
          0.85 0.33 0.10]; % morado
       

% Subplot (a): Estimate from robot 1
subplot(2, 2, 1)
hold on
% Actual trajectories (black dashed)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)


% Estimated trajectories from robot 1 (colored solid)
plot(p1_est1(:,1), p1_est1(:,2), 'Color', colors(1,:), 'LineWidth', 2)


hold off
xlabel('x')
ylabel('y')
title('(a) Estimate 1 from robot 1')
axis equal
grid on
xlim([0 10])
ylim([0 10])

% Subplot (b): Estimate from robot 2
subplot(2, 2, 2)
hold on
% Actual trajectories (black dashed)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)


% Estimated trajectories from robot 2 (colored solid)
plot(p1_est2(:,1), p1_est2(:,2), 'Color', colors(1,:), 'LineWidth', 2)

hold off
xlabel('x')
ylabel('y')
title('(b) Estimate 2 from robot 1')
axis equal
grid on
xlim([0 10])
ylim([0 10])

% Subplot (c): Estimate from robot 3
subplot(2, 2, 3)
hold on
% Actual trajectories (black dashed)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)


% Estimated trajectories from robot 3 (colored solid)
plot(p1_est3(:,1), p1_est3(:,2), 'Color', colors(1,:), 'LineWidth', 2)

hold off
xlabel('x')
ylabel('y')
title('(c) Estimate 3 from robot 1')
axis equal
grid on
xlim([0 10])
ylim([0 10])

% Subplot (f): ESTIMACIÓN PROMEDIADA (Sincronizado de Yu Tang)
subplot(2, 2, 4)
hold on
% Trayectorias reales (negro discontinuo)
plot(p1_actual(:,1), p1_actual(:,2), 'k--', 'LineWidth', 1.5)

% Trayectorias PROMEDIADAS (colores sólidos, más gruesas para destacar)
plot(p1_avg(:,1), p1_avg(:,2), 'Color', colors(1,:), 'LineWidth', 2.5)


hold off
xlabel('x')
ylabel('y')
title('(f) PROMEDIO (Sincronizado)')
axis equal
grid on
xlim([0 10])
ylim([0 10])

%%


figure;
hold on;
plot(t, ERR1, 'Color', colors(1,:), 'LineWidth', 1.2);
plot(t, ERR2, 'Color', colors(2,:), 'LineWidth', 1.2);
plot(t, ERR3, 'Color', colors(3,:), 'LineWidth', 1.2);
plot(t, ERRp, 'k', 'LineWidth', 1.2);

hold off;
xlabel('Tiempo (s)');
ylabel('Error de estimación (norma)');
grid on;
title('Evolución del error de estimación (medición global)');
legend('Observador 1','Observador 2','Observador 3');
%%
figure('Units', 'centimeters', 'Position', [0 0 17 12]);
    plot(t, Z(:,1), 'k', 'LineWidth', 2); hold on
    plot(t, Zob1(:,1), '--', 'Color', colors(1,:), 'LineWidth', 2);
    plot(t, Zob2(:,1), '--', 'Color', colors(2,:), 'LineWidth', 2);
    plot(t, Zob3(:,1), '--', 'Color', colors(3,:), 'LineWidth', 2);
        grid on
    ylabel('Y', 'FontSize', 11, 'Interpreter','latex');
    xlabel('Time [s]','FontSize', 11, 'Interpreter','latex');
    legend('$x$','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','FontSize', 11, 'Interpreter','latex'); % Ajusta los nombres según corresponda
    set(gca, 'FontSize', 11, 'FontName', 'Times', 'LineWidth', 0.5)

    %% 
% ╔═══════════════════════════════════════════════════════════════════════╗
% ║                         CÁLCULO DE ERRORES                            ║
% ╚═══════════════════════════════════════════════════════════════════════╝
error_obs1 = Zob1(:,:) - Z(:,:);  % Error observador 1
error_obs2 = Zob2(:,:) - Z(:,:);  % Error observador 2
error_obs3 = Zob3(:,:) - Z(:,:);  % Error observador 3
error_sinc = Zob_p(:,:) - Z(:,:);  % Error observador sinc

WE11= vecnorm([error_obs1(:,1:2)],  2, 2); %observador 1 del robot 1
WE21= vecnorm([error_obs2(:,1:2)],  2, 2); %observador 2 del robot 1
WE31= vecnorm([error_obs3(:,1:2)],  2, 2); %observador 3 del robot 1
WEsinc= vecnorm([error_sinc(:,1:2)],  2, 2); %observador sinc del robot 1

dt       = ts;
N_robots = 1;
N_obs    = 3;

% Inicializar matrices de métricas
IAE_pos  = zeros(N_obs, N_robots);
ITAE_pos = zeros(N_obs, N_robots);
ISE_pos  = zeros(N_obs, N_robots);

% Llenar para cada observador y robot
% Observador 1
IAE_pos(1,1)  = sum(WE11) * dt;
ITAE_pos(1,1) = sum(t .* WE11') * dt;
ISE_pos(1,1)  = sum(WE11.^2) * dt;

% Observador 2
IAE_pos(2,1)  = sum(WE21) * dt;
ITAE_pos(2,1) = sum(t .* WE21') * dt;
ISE_pos(2,1)  = sum(WE21.^2) * dt;

% Observador 3
IAE_pos(3,1)  = sum(WE31) * dt;
ITAE_pos(3,1) = sum(t .* WE31') * dt;
ISE_pos(3,1)  = sum(WE31.^2) * dt;

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
    round(promedio_ITAE_obs, 3), ...
    round(promedio_IAE_obs, 3), ...
    round(promedio_ISE_obs, 3), ...
    'VariableNames', {'ITAE Average', 'IAE Average', 'ISE Average'}, ...
    'RowNames', {'x1_hat', 'x2_hat', 'x3_hat'}); 
disp('=== PROMEDIO DE MÉTRICAS ===')
disp(tabla_metricas)