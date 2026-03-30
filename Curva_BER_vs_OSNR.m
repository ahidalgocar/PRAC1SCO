addpath('D:\Program Files\MATLAB\optilux\Optilux_files') 

clear all; clc; close all;

% Parámetros de la simulación
Nsymb = 512;   
Nsym = 512;
Nt = 32;      
Nch = 1;       

% Parámetros de la señal
phi = 0.2*pi; 
exratio = 20; 
extratio=20;
spac = 0.4;   
symbrate = 10.7; 
duty = 1;
roll = 0.2;

% Parámetros del enlace
Din = 0;       
Nspan = 83;   

% Variables para guardar resultados
pb_1530 = [];
pb_1625 = [];
osnr_ref = [];
OSNR_B2B = [];

% Bucle para ambas longitudes de onda
for lam = [1530 1625]

    % Parámetros dependientes de la longitud de onda
    if lam == 1530
        tx.alphadB = 0.18;     
        tx.disp = 16.1;        
        comp.disp = -139.5;    
    else
        tx.alphadB = 0.2;
        tx.disp = 21.4;
        comp.disp = -186;
    end

    % Fibra de transmisión
    tx.length = 72000;   
    tx.aeff = 80;        
    tx.n2 = 2.7e-20;     
    tx.lambda = lam;     
    tx.slope = 0;
    tx.dphimax = 3E-3;
    tx.dzmax = 2E4;

    % Fibra compensadora
    comp.alphadB = 0.4;
    comp.aeff = 20;
    comp.n2 = 2.7e-20;
    comp.lambda = lam;
    comp.slope = 0;
    comp.dphimax = 3E-3;
    comp.dzmax = 2E4;

    % Receptor
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

    % Conversión
    gam = 2*pi*tx.n2./(lam.*tx.aeff)*1e21;

    % Longitud DCF
    comp.length = (Din - tx.disp*tx.length/1e3)/comp.disp*1e3;

    % Ganancia EDFA
    Gerbio = tx.alphadB*tx.length/1e3 + comp.alphadB*comp.length/1e3;
    ampli.f = -1;

    % Generación señal
    reset_all(Nsymb,Nt,Nch);

    Pavg = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan);

    E = lasersource(Pavg, lam, spac);

    for ii=1:Nch
        pat(:,ii)=olpattern('debruijn',ii*2);
        elec(:,ii)=electricsource(pat(:,ii),'ook',symbrate,'cosroll',duty,roll);
        Eopt(:,ii) = mz_modulator(E(:,ii),elec(:,ii),struct('exratio',exratio));
    end

    create_field('unique',Eopt,[],struct('power','average'));

    % BACK-TO-BACK
    x.b2b = 'b2b';
    [pb_b2b, osnr_b2b] = ber_kl(1,x,pat);
    x = rmfield(x,'b2b');

    % PROPAGACIÓN 
    for k=1:Nspan
        fiber(tx,'g-sx')
        fiber(comp,'g-sx')
        ampliflat(Gerbio,'gain', ampli)
    end

    % BER vs OSNR
    [pb,osnr]=ber_kl(1,x,pat);

    % Guardar resultados
    if lam == 1530
        pb_1530 = pb;
        osnr_ref = x.osnr;
        OSNR_B2B_1530 = osnr_b2b;
    else
        pb_1625 = pb;
        OSNR_B2B_1625 = osnr_b2b;
    end

    % Figura individual
    figure
    semilogy(x.osnr,pb);
    grid on;
    xlabel('OSNR [dB/0.1nm]')
    ylabel('BER')
    title(['BER vs OSNR - ', num2str(lam), ' nm'])

end

% FIGURA COMPARATIVA 
figure
semilogy(osnr_ref, pb_1530, 'LineWidth', 2)
hold on
semilogy(osnr_ref, pb_1625, 'LineWidth', 2)

grid on
xlabel('OSNR [dB/0.1nm]')
ylabel('BER')
title('BER vs OSNR (comparación)')
legend('1530 nm','1625 nm')

ylim([1e-5 1e-1]) 

% CORTES BER = 1e-3 
BER_obj = 1e-3;

OSNR_1530 = interp1(log10(pb_1530), osnr_ref, log10(BER_obj));
OSNR_1625 = interp1(log10(pb_1625), osnr_ref, log10(BER_obj));

%  PENALIZACIÓN 
pen_1530 = OSNR_1530 - OSNR_B2B_1530;
pen_1625 = OSNR_1625 - OSNR_B2B_1625;

% MOSTRAR RESULTADOS 
fprintf('\n--- RESULTADOS ---\n')
fprintf('B2B 1530 nm: %.2f dB\n', OSNR_B2B_1530)
fprintf('B2B 1625 nm: %.2f dB\n', OSNR_B2B_1625)

fprintf('\n--- Corte BER = 1e-3 ---\n')
fprintf('1530 nm: %.2f dB\n', OSNR_1530)
fprintf('1625 nm: %.2f dB\n', OSNR_1625)

fprintf('\n--- Penalización ---\n')
fprintf('1530 nm: %.2f dB\n', pen_1530)
fprintf('1625 nm: %.2f dB\n', pen_1625)

% DIBUJO EN FIGURA
yline(BER_obj,'k--','BER = 10^{-3}','LineWidth',1.5);

plot(OSNR_1530, BER_obj, 'ro', 'MarkerSize',8,'LineWidth',2)
plot(OSNR_1625, BER_obj, 'bo', 'MarkerSize',8,'LineWidth',2)

text(OSNR_1530, BER_obj*1.5, sprintf('%.2f dB',OSNR_1530), 'Color','r')
text(OSNR_1625, BER_obj*1.5, sprintf('%.2f dB',OSNR_1625), 'Color','b')