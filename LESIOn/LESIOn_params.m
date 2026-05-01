%%%% Frankenmodel Parameters %%%%

function p = LESIOn_params()

% Miller Parameters
    p.beta1 = 1e-6;
    p.omega = 1e-5;
    p.mum = 1e4;
    p.muk = 1e1; % old: 1e4, new: 3.5e3
    p.etaE = 1;
    p.etaK = 0.02;
    p.etaM = 0.2;
    p.deltaM0 = 0.02;
    p.deltaM1 = 0.02;
    p.deltaM2 = 0.02;
    p.deltaK0 = 0.02;
    p.deltaKa = 0.02;
    p.deltaE0 = 1;
    p.deltaEa = 0.14;
    p.deltaEf = 0.14;
    p.sigma = 1e-8 ; % old: 1e-5
    p.mc = 1e6;
    p.kc = 1e6;
    p.ec = 1e9;
    p.cm1 = 1e5;
    p.cm2 = 100;
    p.cka = 1e5; 
    p.beta2 = 1e-3;
    p.beta12 = 5e-5;
    p.beta21 = 5e-5;
    p.thetaM = 5e-8;
    p.thetaK = 0.7; % old: 0.2, new: 0.7;
    p.gamma = 0.8;
    p.rho0 = 0.1;
    p.rhoF = 0.1;
   % mue = 287.6;
    p.mue = 1000;
 % Hormonal modifications
    % NK cells
    p.phi1 = 20; % old: phi1 = 10 ;
   % n1 = 1; % old: n1 = 1;
    p.H1 = 1;
    % endometrial cells
    p.phi2 = 100;
    p.n2 = 1;
    p.H2 = 10;
    % receptor binding
    %kalphaM = 1.48e-2 * 60 * 60 * 24; % szatkowski 2001, 1/s
    p.kalphaM = 1.2e-3 * 60 * 60 * 24 ; % benchchem
    p.kalphaP = 1.3e6 * 60 * 60 * 24; % benchchem, 1/M*s
    p.kbetaM = 3.0e-3 * 60 * 60 * 24; % benchchem, 1/s
    p.kbetaP = 1.0e6 * 60 * 60 * 24; % benchchem, 1/M*s
    p.calphaM = 3.0e-4 * 60 * 60 * 24; % calculated from half life (koff = ln(2) / t1/2) - tamrazi 2002
    p.calphaP = 1.5e4 * 60 * 60 * 24; % jisa 2003, 1/M*s
    p.cbetaM = 2.2e-4 * 60 * 60 * 24; % jisa 2003 - est. Kd / kon = koff
    p.cbetaP = 5.7e3 * 60 * 60 * 24; % jisa 2003, 1/M*s
    p.phi3 = 100;
    p.n3 = 10;
    p.H3 = 1;
    p.phi4 = 100;
    p.n4 = 1;
    p.H4 = 10;
    p.phi5 = 10;
    p.n5 = 8;
    p.H5 = 0.25;
  % Graham Parameters
    % pituitary parameters
  p.vf = 3219.9 ;
  p.KiFI = 149.76 ;
  p.kF = 3.0212 ;
  p.cfp = 65.229 ;
  p.cfe = 0.0024047 ;
  p.cfi = 3.0188 ;
  p.v0L = 308.35 ;
  p.v1L = 44700 ;
  p.KmL = 226.37 ;
  p.KiLP = 3.2279 ;
  p.kL = 0.67146 ;
  p.clp = 0.015844 ;
  p.cle = 0.00068867 ;

  % ovarian parameters
  p.f1 = 1.0958 ;
  p.f2 = 46.225 ;
  p.h1 = 146.31 ;
  p.h2 = 798.39 ;
  p.w = 0.23497 ;
  p.l = 0.64178 ;
  p.shat = 2.6338 ;
  p.deltaS = 0.38256 ;
  p.eta = 0.81426 ;
  p.k2 = 8.276 ;
  p.hs = 11.691 ;
  p.tg1 = 6.3594 ;
  p.ee0 = 9.6377 ;
  p.p = 0.22851 ;

  % fixed parameters
  p.deltaF = 8.21 ;
  p.V = 2.5 ;
  p.deltaL = 14;
  p.deltaE = 1.1;
  p.deltaP = 0.5 ;

  % from selgrade
  p.n = 8;
  
  % local E2
  p.deltaLE2 = 0.03;
  p.pk = 0.005;
end
