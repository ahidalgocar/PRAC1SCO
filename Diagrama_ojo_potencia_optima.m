addpath('D:\Program Files\MATLAB\optilux\Optilux_files')

clear all; clc; close all;

%% Parámetros comunes
Nsymb = 512;
Nt = 32;
Nch = 1;

symbrate = 10.7;
exratio = 20;
Din = 0;
Nspan = 83;

%% ============================
%% DEFINICIÓN DE LOS 6 CASOS
%% ============================

% 1530 nm
casos_1530 = [ ...
    -5   0.15*pi;   % mínimo
    -1   0.20*pi;   % óptimo
     2   0.30*pi];  % máximo

% 1625 nm
casos_1625 = [ ...
    -5   0.13*pi;   % mínimo
    -2   0.18*pi;   % óptimo
     2   0.28*pi];  % máximo

%% ============================
%% FUNCIÓN PRINCIPAL
%% ============================

for lam = [1530 1625]

    if lam == 1530
        casos = casos_1530;
        tx.alphadB = 0.18;
        tx.disp = 16.1;
        comp.disp = -139.5;
    else
        casos = casos_1625;
        tx.alphadB = 0.2;
        tx.disp = 21.4;
        comp.disp = -186;
    end

    %% Fibra transmisión
    tx.length = 72000;
    tx.lambda = lam;
    tx.aeff = 80;
    tx.n2 = 2.7e-20;
    tx.slope = 0;
    tx.dphimax = 3E-3;
    tx.dzmax = 2E4;

    %% Fibra compensadora
    comp.alphadB = 0.4;
    comp.lambda = lam;
    comp.aeff = 20;
    comp.n2 = 2.7e-20;
    comp.slope = 0;
    comp.dphimax = 3E-3;
    comp.dzmax = 2E4;

    comp.length = (Din - tx.disp*tx.length/1e3)/comp.disp*1e3;

    %% Amplificador
    Gerbio = tx.alphadB*tx.length/1e3 + comp.alphadB*comp.length/1e3;
    ampli.f = -1;

    %% Parámetro no lineal
    gam = 2*pi*tx.n2/(lam*tx.aeff)*1e21;

    %% ============================
    %% BUCLE SOBRE LOS 3 CASOS
    %% ============================

    for i = 1:size(casos,1)

        P_dBm = casos(i,1);
        phi   = casos(i,2);

        %% Reset sistema
        reset_all(Nsymb,Nt,Nch);

        %% Potencia (CLAVE)
        Pavg = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan);

        %% Generación señal
        E = lasersource(Pavg, lam, 0.4);

        pat = olpattern('debruijn',2);
        elec = electricsource(pat,'ook',symbrate,'cosroll',1,0.2);
        Eopt = mz_modulator(E,elec,struct('exratio',exratio));

        create_field('unique',Eopt,[],struct('power','average'));

        %% Propagación
        for k = 1:Nspan
            fiber(tx,'g-sx')
            fiber(comp,'g-sx')
            ampliflat(Gerbio,'gain', ampli)
        end

        %% ============================
        %% DIAGRAMA DE OJO
        %% ============================

        figure;
        grid on; hold on;

        x.oftype = 'gauss';
        x.obw = 1.8;
        x.eftype = 'bessel5';
        x.ebw = 0.65;
        x.rec = 'ook';
        x.plot = 'ploteye';
        x.color = 'b-';

        x.ber = 1e-3;
        x.eta = 1.4;
        x.mu = 3.5;
        x.osnr = 15+(-7:10);
        x.poln = 2;
        x.saddle = 'y';

        % SOLO PROPAGACIÓN (sin b2b para evitar errores)
        [pb,osnr] = ber_kl(1,x,pat);

        title(['Diagrama de ojo - ', num2str(lam), ...
               ' nm @ ', num2str(P_dBm), ' dBm']);

    end
end