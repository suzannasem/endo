%%% File containing parameters for Graham Model %%%
function [vf, KiFI, kF, cfp, cfe, cfi, v0L, v1L, KmL, KiLP, ...
    kL, clp, cle, f1, f2, h1, h2, w, l, shat, deltaS, eta, k2, ...
    hs, tg1, e0, p, deltaF, V, deltaL, deltaE, deltaP, ...
    n] = graham_params
  
% pituitary parameters
  vf = 3219.9 ;
  KiFI = 149.76 ;
  kF = 3.0212 ;
  cfp = 65.229 ;
  cfe = 0.0024047 ;
  cfi = 3.0188 ;
  v0L = 308.35 ;
  v1L = 44700 ;
  KmL = 226.37 ;
  KiLP = 3.2279 ;
  kL = 0.67146 ;
  clp = 0.015844 ;
  cle = 0.00068867 ;

  % ovarian parameters
  f1 = 1.0958 ;
  f2 = 46.225 ;
  h1 = 146.31 ;
  h2 = 798.39 ;
  w = 0.23497 ;
  l = 0.64178 ;
  shat = 2.6338 ;
  deltaS = 0.38256 ;
  eta = 0.81426 ;
  k2 = 8.276 ;
  hs = 11.691 ;
  tg1 = 6.3594 ;
  e0 = 9.6377 ;
  p = 0.22851 ;

  % fixed parameters
  deltaF = 8.21 ;
  V = 2.5 ;
  deltaL = 14;
  deltaE = 1.1;
  deltaP = 0.5 ;

  % from selgrade
  n = 8;

end
