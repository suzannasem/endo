function [mum,muk,etaE,etaK,etaM,deltaM0,deltaM1,deltaM2,deltaK0,deltaKa,deltaE0,deltaEa, deltaEf,sigma,mc, ...
              kc,ec,cka,cm1,cm2,beta2,beta12,beta21,thetaM,thetaK,gamma, ...
                 rhoF, rho0, mue,...
                 phi1, n1, H1, phi2, n2, H2, ...
                 kalphaM, kalphaP, kbetaM, kbetaP, calphaM, calphaP, cbetaM, cbetaP, ... % K1, K2, ...
                 phi3, n3, H3, phi4, n4, H4, ...
                 phi5, n5, H5] = Constant_E2_params

    mum = 1e4;
    muk = 1e1; % old: 1e4, new: 3.5e3
    etaE = 1;
    etaK = 0.02;
    etaM = 0.2;
    deltaM0 = 0.02;
    deltaM1 = 0.02;
    deltaM2 = 0.02;
    deltaK0 = 0.02;
    deltaKa = 0.02;
    deltaE0 = 1;
    deltaEa = 0.14;
    deltaEf = 0.14;
    %sigma = 1e-5;
    sigma = 1e-8 ;
    mc = 1e6;
    kc = 1e6;
    ec = 1e9;
    cm1 = 1e5;
    cm2 = 100;
    cka = 1e5; % old: 1e6, new: 1e5
    beta2 = 1e-3;
    beta12 = 5e-5;
    beta21 = 5e-5;
    thetaM = 5e-8;
    thetaK = 0.7; % old: 0.2, new: 0.7;
    gamma = 0.8;
    rho0 = 0.1;
    rhoF = 0.1;
    mue = 1260;
   % Hormonal modifications
    % NK cells
    phi1 = 10;
    n1 = 1;
    H1 = 1;
    % endometrial cells
    phi2 = 100;
    n2 = 1;
    H2 = 10;
    % macrophages
    %kalphaM = 1.48e-2 * 60 * 60 * 24; % szatkowski 2001, 1/s
    kalphaM = 1.2e-3 * 60 * 60 * 24 ; % benchchem
    kalphaP = 1.3e6 * 60 * 60 * 24; % benchchem, 1/M*s
    kbetaM = 3.0e-3 * 60 * 60 * 24; % benchchem, 1/s
   % kalphaM = 1.2e-3 * 60 * 60 * 24;
   % kalphaP = 1.3e7 * 60 * 60 * 24;
   % kbetaM = 3.0e3 * 60 * 60 * 24;
   % kbetaP = 1.0e-6 * 60 * 60 * 24;
    kbetaP = 1.0e6 * 60 * 60 * 24; % benchchem, 1/M*s
    calphaM = 3.0e-4 * 60 * 60 * 24; % calculated from half life (koff = ln(2) / t1/2) - tamrazi 2002
    calphaP = 1.5e4 * 60 * 60 * 24; % jisa 2003, 1/M*s
    cbetaM = 2.2e-4 * 60 * 60 * 24; % jisa 2003 - est. Kd / kon = koff
    cbetaP = 5.7e3 * 60 * 60 * 24; % jisa 2003, 1/M*s
  %  K1 = 5*10^4; % guess
  %  K2 = 3*10^3;
    phi3 = 100;
    n3 = 10;
    H3 = 1;
    phi4 = 100;
    n4 = 1;
    H4 = 10;
    phi5 = 10;
    n5 = 8;
    H5 = 0.25;
end

