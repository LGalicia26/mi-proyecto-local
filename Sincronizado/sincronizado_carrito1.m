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

z = [2.5;  0.5; 0; 0];   
 
% Initial conditions for observers - diferentes de los reales
zob1 = z + 0.4*randn(4,1);
zob2 = z + 0.4*randn(4,1);
zob3 = z + 0.4*randn(4,1);


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

%% Tiempos de muestreo
H1 = 0.04; 
H2 = H1*2;  
H3 = H1*3;  

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
    u = -K * [zeros(2,1); zob1(3:4)  - vel_deseada];
    dz = A*z + B*u;

     %% Muestreo 1
        rui1=0.001;
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
        rui2=0.002;
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
        rui3=0.003;
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
    
    %% Estimación promedio de las estimaciones. 

    % Estimación promedio de las estimaciones de posición y velocidad  
        zob_p = (zob1+zob2+zob3)/3;
        Zob_p = [Zob_p;zob_p'];

    % Error de la estimación promedio de posición y velocidad
        erp    = norm(z - zob_p);
        ERRp   = [ERRp; erp'];

        j=j+1;
end
% ╔═══════════════════════════════════════════════════════════════════════╗
% ║                         CÁLCULO DE ERRORES                            ║
% ╚═══════════════════════════════════════════════════════════════════════╝
error_obs1 = Zob1(:,:) - Z(:,:);   % Error observador 1
error_obs2 = Zob2(:,:) - Z(:,:);   % Error observador 2
error_obs3 = Zob3(:,:) - Z(:,:);   % Error observador 3
error_sinc = Zob_p(:,:) - Z(:,:);  % Error observador sinc

WE11   = vecnorm(error_obs1(:,1:2),  2, 2); % observador 1 del robot 1
WE21   = vecnorm(error_obs2(:,1:2),  2, 2); % observador 2 del robot 1
WE31   = vecnorm(error_obs3(:,1:2),  2, 2); % observador 3 del robot 1
WEsinc = vecnorm(error_sinc(:,1:2),  2, 2); % observador sinc del robot 1

dt       = ts;
N_robots = 1;
N_obs    = 4; % Cambiado a 4 para incluir el sincronizado

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

% Observador Sincronizado (Promedio)
IAE_pos(4,1)  = sum(WEsinc) * dt;
ITAE_pos(4,1) = sum(t .* WEsinc') * dt;
ISE_pos(4,1)  = sum(WEsinc.^2) * dt;

% Promedios por observador
promedio_IAE_obs  = mean(IAE_pos, 2); 
promedio_ITAE_obs = mean(ITAE_pos, 2); 
promedio_ISE_obs  = mean(ISE_pos, 2); 

% (Opcional si usas más robots) Promedios por robot
promedio_IAE_robot  = mean(IAE_pos, 1); 
promedio_ITAE_robot = mean(ITAE_pos, 1); 
promedio_ISE_robot  = mean(ISE_pos, 1); 

%%  TABLA CON ERRORES
tabla_metricas = table(...
    round(promedio_ITAE_obs, 3), ...
    round(promedio_IAE_obs, 3), ...
    round(promedio_ISE_obs, 3), ...
    'VariableNames', {'ITAE Average', 'IAE Average', 'ISE Average'}, ...
    'RowNames', {'x1_hat', 'x2_hat', 'x3_hat', 'x_sinc'}); % Agregada la fila x_sinc

disp('=== PROMEDIO DE MÉTRICAS ===')
disp(tabla_metricas)

% Define colors
colors = [0.71 0.43 0.47; % rosa
          0.46 0.67 0.19; % verde
          0.85 0.33 0.10]; % morado
       
%% GRÁFICA DE ERROR DE ESTIMACIÓN (CON VENTANA DE ZOOM)

% 1. Definir el intervalo de tiempo para el acercamiento (Zoom)
% Elegimos un intervalo en estado estacionario (ej. de 10 a 14 segundos)
t_zoom_err = [10 14]; 
idx_err = find(t >= t_zoom_err(1) & t <= t_zoom_err(2));

% Encontrar el valor máximo del error en esa ventana para dibujar el rectángulo
% (El mínimo será cero porque es una norma)
y_max_err = max([max(ERR1(idx_err)), max(ERR2(idx_err)), max(ERR3(idx_err)), max(ERRp(idx_err))]);
y_min_err = 0; 

% Crear la figura principal
figure('Name', 'Error de estimación', 'Position', [50, 50, 600, 400]);

% --- Gráfica Principal ---
hold on;
plot(t, ERR1, 'Color', colors(1,:), 'LineWidth', 1.2);
plot(t, ERR2, 'Color', colors(2,:), 'LineWidth', 1.2);
plot(t, ERR3, 'Color', colors(3,:), 'LineWidth', 1.2);
plot(t, ERRp, 'k', 'LineWidth', 1.5); % Sincronizado un poco más grueso

% Dibujar rectángulo rojo indicando la zona del zoom
% Le sumamos un pequeño margen al y_max_err para que el cuadro no quede tan apretado
margen = y_max_err * 0.2; 
rectangle('Position', [t_zoom_err(1), y_min_err, t_zoom_err(2)-t_zoom_err(1), y_max_err+margen], ...
          'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '-.');
hold off;

% Formato de la gráfica principal
xlabel('Tiempo (s)', 'FontSize', 12);
ylabel('Error de estimación (norma)', 'FontSize', 12);
grid on;
title('Evolución del error de estimación (medición global)', 'FontSize', 13);
legend('Observador 1', 'Observador 2', 'Observador 3', 'Sincronizado', ...
       'Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 11, 'FontName', 'Times');

% --- Ventana de Zoom (Inset) ---
% [izquierda, abajo, ancho, alto] - Ajustado para no tapar el pico inicial
axes('Position', [0.45 0.4 0.4 0.35]); 
box on; hold on;
plot(t(idx_err), ERR1(idx_err), 'Color', colors(1,:), 'LineWidth', 1.2);
plot(t(idx_err), ERR2(idx_err), 'Color', colors(2,:), 'LineWidth', 1.2);
plot(t(idx_err), ERR3(idx_err), 'Color', colors(3,:), 'LineWidth', 1.2);
plot(t(idx_err), ERRp(idx_err), 'k', 'LineWidth', 1.5);
hold off;

% Formato de la ventana de zoom
grid on;
title('Detalle en estado estacionario', 'FontSize', 10);
xlim(t_zoom_err);
% Fijar límites en Y para la ventana usando el máximo calculado
ylim([0, y_max_err + margen]); 
set(gca, 'FontSize', 9, 'FontName', 'Times');

%% GRÁFICAS COMPARATIVAS: Ruido vs Estimaciones (CON VENTANAS DE ZOOM)

% 1. Definir el intervalo de tiempo para el acercamiento (Zoom)
% Reducimos el zoom a un segundo (de 10 a 11) para que el ruido se vea mucho más claro
t_zoom = [10 11]; 
idx_zoom = find(t >= t_zoom(1) & t <= t_zoom(2));

% Calcular los límites Y del rectángulo del zoom para dibujarlo correctamente
y_min_zoom = min(Y3(idx_zoom, 1));
y_max_zoom = max(Y3(idx_zoom, 1));

% Crear figura grande para acomodar los insets
figure('Name', 'Efecto del Ruido: Principal con Zoom', 'Position', [50, 50, 600, 400]);

% ╔═══════════════════════════════════════════════════════════════════════╗
% ║       GRÁFICA 1: Ruido con Estimaciones Individuales (1, 2 y 3)       ║
% ╚═══════════════════════════════════════════════════════════════════════╝
% --- Gráfica Principal ---
subplot(2, 1, 1)
hold on;
% Medición ruidosa (Gris fuerte y sólida)
plot(t, Y3(:, 1), '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
plot(t, Z(:, 1), 'k', 'LineWidth', 2);
plot(t, Zob1(:, 1), '--', 'Color', colors(1,:), 'LineWidth', 1.5);
plot(t, Zob2(:, 1), '--', 'Color', colors(2,:), 'LineWidth', 1.5);
plot(t, Zob3(:, 1), '--', 'Color', colors(3,:), 'LineWidth', 1.5);

% Dibujar un rectángulo rojo indicando la zona del zoom
rectangle('Position', [t_zoom(1), y_min_zoom, t_zoom(2)-t_zoom(1), y_max_zoom-y_min_zoom], ...
          'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '-.');
hold off;
grid on;
ylabel('Posición X', 'FontSize', 12, 'Interpreter', 'latex');
title('Ruido vs Estimaciones Individuales (Trayectoria Completa)', 'FontSize', 13);
legend('Ruido ($y_3$)', 'Real ($x$)', '$\hat{x}_1$', '$\hat{x}_2$', '$\hat{x}_3$', ...
       'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 11);
set(gca, 'FontSize', 11, 'FontName', 'Times');

% --- Ventana de Zoom (Inset) 1 ---
% [izquierda, abajo, ancho, alto] en coordenadas normalizadas de la figura
axes('Position', [0.18 0.72 0.25 0.15]); 
box on; hold on;
plot(t(idx_zoom), Y3(idx_zoom, 1), '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
plot(t(idx_zoom), Z(idx_zoom, 1), 'k', 'LineWidth', 2);
plot(t(idx_zoom), Zob1(idx_zoom, 1), '--', 'Color', colors(1,:), 'LineWidth', 1.5);
plot(t(idx_zoom), Zob2(idx_zoom, 1), '--', 'Color', colors(2,:), 'LineWidth', 1.5);
plot(t(idx_zoom), Zob3(idx_zoom, 1), '--', 'Color', colors(3,:), 'LineWidth', 1.5);
hold off;
grid on;
xlim(t_zoom);
set(gca, 'FontSize', 9, 'FontName', 'Times');

% ╔═══════════════════════════════════════════════════════════════════════╗
% ║            GRÁFICA 2: Ruido con Estimación Sincronizada               ║
% ╚═══════════════════════════════════════════════════════════════════════╝
% --- Gráfica Principal ---
subplot(2, 1, 2)
hold on;
plot(t, Y3(:, 1), '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
plot(t, Z(:, 1), 'k', 'LineWidth', 2);
plot(t, Zob_p(:, 1), '-', 'Color', 'b', 'LineWidth', 2);

% Dibujar rectángulo rojo
rectangle('Position', [t_zoom(1), y_min_zoom, t_zoom(2)-t_zoom(1), y_max_zoom-y_min_zoom], ...
          'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '-.');
hold off;
grid on;
xlabel('Tiempo [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Posición X', 'FontSize', 12, 'Interpreter', 'latex');
title('Ruido vs Estimación Sincronizada (Trayectoria Completa)', 'FontSize', 13);
legend('Ruido ($y_3$)', 'Real ($x$)', '$\hat{x}_{sinc}$', ...
       'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 11);
set(gca, 'FontSize', 11, 'FontName', 'Times');

% --- Ventana de Zoom (Inset) 2 ---
axes('Position', [0.18 0.25 0.25 0.15]);
box on; hold on;
plot(t(idx_zoom), Y3(idx_zoom, 1), '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
plot(t(idx_zoom), Z(idx_zoom, 1), 'k', 'LineWidth', 2);
plot(t(idx_zoom), Zob_p(idx_zoom, 1), '-', 'Color', 'b', 'LineWidth', 2);
hold off;
grid on;
xlim(t_zoom);
set(gca, 'FontSize', 9, 'FontName', 'Times');