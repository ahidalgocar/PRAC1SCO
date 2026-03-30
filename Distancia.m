clc; clear; close all;

%% Parámetros
P0 = 2e-3;                 % Potencia [W]
Fn = 10^(5.5/10);          % Factor de ruido
h = 6.626e-34;             % Constante de Planck
B0 = 12.5e9;               % Ancho de banda [Hz]
c = 3e8;                   % Velocidad luz

%% Longitudes de onda
lambda1 = 1530e-9;
lambda2 = 1625e-9;

nu1 = c/lambda1;
nu2 = c/lambda2;

%% Parámetros de fibra
alpha_SMF_1530 = 0.18; % dB/km
alpha_SMF_1625 = 0.2;  % dB/km
alpha_DCF = 0.4;       % dB/km

%% Diseño
L_SMF = 72; % km
L_DCF = 8;  % km
L_span = L_SMF + L_DCF;

N = ceil(6605 / L_span);

%% Pérdidas por tramo (dB)
L_1530_dB = alpha_SMF_1530*L_SMF + alpha_DCF*L_DCF;
L_1625_dB = alpha_SMF_1625*L_SMF + alpha_DCF*L_DCF;

%% Conversión a lineal
L_1530 = 10^(L_1530_dB/10);
L_1625 = 10^(L_1625_dB/10);

%% Cálculo OSNR
OSNR_1530 = P0 / (N * Fn * h * nu1 * L_1530 * B0);
OSNR_1625 = P0 / (N * Fn * h * nu2 * L_1625 * B0);

%% Pasar a dB
OSNR_1530_dB = 10*log10(OSNR_1530);
OSNR_1625_dB = 10*log10(OSNR_1625);

%% Mostrar resultados
fprintf('Número de tramos: %d\n', N);

fprintf('\n--- 1530 nm ---\n');
fprintf('Pérdidas tramo: %.2f dB\n', L_1530_dB);
fprintf('OSNR: %.2f (%.2f dB)\n', OSNR_1530, OSNR_1530_dB);

fprintf('\n--- 1625 nm ---\n');
fprintf('Pérdidas tramo: %.2f dB\n', L_1625_dB);
fprintf('OSNR: %.2f (%.2f dB)\n', OSNR_1625, OSNR_1625_dB);