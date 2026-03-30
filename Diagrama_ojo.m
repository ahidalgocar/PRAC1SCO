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

%% Bucle para ambas longitudes de onda
for lam = [1530 1625]

    %% Fibra SMF
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

    %% Señal
    reset_all(Nsymb,Nt,Nch);

    phi = 0.2*pi;
    gam = 2*pi*tx.n2/(lam*tx.aeff)*1e21;

    Pavg = phi2pow(phi,tx.length,tx.alphadB,gam,0,Nspan);

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

    %% DIAGRAMA DE OJO
    figure;
    grid on; hold on;

    x.oftype = 'gauss';
    x.obw = 1.8;
    x.eftype = 'bessel5';
    x.ebw = 0.65;
    x.rec = 'ook';
    x.plot = 'ploteye';
    x.color = 'r-';

    x.ber = 1e-3;
    x.eta = 1.4;
    x.mu = 3.5;
    x.osnr = 15+(-7:10);
    x.poln = 2;
    x.saddle = 'y';

    % Back-to-back
    x.b2b = 'b2b';
    [~,~] = ber_kl(1,x,pat);

    % Propagación
    x = rmfield(x,'b2b');
    x.color = 'b-';
    [~,~] = ber_kl(1,x,pat);

    title(['Diagrama de ojo - ', num2str(lam), ' nm']);

end