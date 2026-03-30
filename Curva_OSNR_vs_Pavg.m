addpath('D:\Program Files\MATLAB\optilux\Optilux_files')

clear all; clc; close all;

% Parámetros del campo
Nsymb = 512;
Nsym = 512;
Nt = 32;
Nch = 1;

% Parámetros de la señal
phi0 = 0.2*pi;
exratio = 20;
extratio=20;
spac = 0.4;
symbrate = 10.7;
duty = 1;
roll = 0.2;

% Parámetros del enlace
Din = 0;
Nspan = 83;

% Parámetros del receptor
x.oftype = 'gauss';
x.obw = 1.8;
x.eftype = 'bessel5';
x.ebw = 0.65;
x.rec = 'ook';

x.ber = 1e-3;
x.eta = 1.4;
x.mu = 3.5;
x.osnr = 15+(-7:10);
x.poln = 2;
x.saddle = 'y';

% Vector de potencia
phi_vec = 0.1*pi : 0.05*pi : 0.5*pi;

OSNR_1530 = zeros(size(phi_vec));
OSNR_1625 = zeros(size(phi_vec));

for lambda_case = [1530 1625]

    lam = lambda_case;

    if lam == 1530
        tx.alphadB = 0.18;
        tx.disp = 16.1;
        comp.disp = -139.5;
    else
        tx.alphadB = 0.2;
        tx.disp = 21.4;
        comp.disp = -186;
    end

    tx.length = 72000;
    tx.aeff = 80;
    tx.n2 = 2.7e-20;
    tx.lambda = lam;
    tx.slope = 0;
    tx.dphimax = 3E-3;
    tx.dzmax = 2E4;

    comp.alphadB = 0.4;
    comp.aeff = 20;
    comp.n2 = 2.7e-20;
    comp.lambda = lam;
    comp.slope = 0;
    comp.dphimax = 3E-3;
    comp.dzmax = 2E4;

    gam = 2*pi*tx.n2./(lam.*tx.aeff)*1e21;

    comp.length = (Din - tx.disp*tx.length/1e3)/comp.disp*1e3;

    Gerbio = tx.alphadB*tx.length/1e3 + comp.alphadB*comp.length/1e3;
    ampli.f = -1;

    for i = 1:length(phi_vec)

        phi = phi_vec(i);

        reset_all(Nsymb,Nt,Nch);

        Pavg = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan);

        E = lasersource(Pavg, lam, spac);

        for ii=1:Nch
            pat(:,ii)=olpattern('debruijn',ii*2);
            elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
            Eopt(:,ii) = mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
        end

        create_field('unique',Eopt,[],struct('power','average'));

        for k=1:Nspan
            fiber(tx,'g-sx')
            fiber(comp,'g-sx')
            ampliflat(Gerbio,'gain', ampli)
        end

        [pb,osnr]=ber_kl(1,x,pat);

        if lam == 1530
            OSNR_1530(i) = osnr;
        else
            OSNR_1625(i) = osnr;
        end

    end
end

% Conversión a dBm
P_dBm = 10*log10(phi_vec);

% Gráfica
figure
plot(P_dBm, OSNR_1530,'-o','LineWidth',2)
hold on
plot(P_dBm, OSNR_1625,'-o','LineWidth',2)

grid on
xlabel('Pavg [dBm]')
ylabel('OSNR [dB/0.1 nm]')
title('OSNR vs Potencia (BER = 10^{-3})')
legend('1530 nm','1625 nm')

% =========================
% TABLA DE RESULTADOS
% =========================

Tabla_Resultados = table(P_dBm', OSNR_1530', OSNR_1625', ...
    'VariableNames', {'Potencia_dBm','OSNR_1530_dB','OSNR_1625_dB'});

disp(' ')
disp('--- TABLA DE RESULTADOS ---')
disp(Tabla_Resultados)

% =========================
% POTENCIA ÓPTIMA
% =========================

[OSNR_min_1530, idx1] = min(OSNR_1530);
[OSNR_min_1625, idx2] = min(OSNR_1625);

fprintf('\n--- POTENCIA ÓPTIMA ---\n')
fprintf('1530 nm: %.2f dBm (OSNR = %.2f dB)\n', P_dBm(idx1), OSNR_min_1530)
fprintf('1625 nm: %.2f dBm (OSNR = %.2f dB)\n', P_dBm(idx2), OSNR_min_1625)