%%%% Frankenmodel: Graham (2023) + Miller (2025) + Semaan (2026) %%%%

function out = Franken_fn(t,Y, beta1, omega, K1, K2, psi)
%% parameters %%

[mum,muk,etaE,etaK,etaM,deltaM0,deltaM1,deltaM2,deltaK0,deltaKa,deltaE0,deltaEa, deltaEf,sigma,mc, ...
              kc,ec,cka,cm1,cm2,beta2,beta12,beta21,thetaM,thetaK,gamma, ...
                 rhoF, rho0, mue, ...
                 phi1, H1, phi2, n2, H2, ...
                 kalphaM, kalphaP, kbetaM, kbetaP, calphaM, calphaP, cbetaM, cbetaP, ... % K1, K2, ...
                 phi3, n3, H3, phi4, n4, H4, ...
                 phi5, n5, H5, ...
                 vf, KiFI, kF, cfp, cfe, cfi, v0L, v1L, KmL, KiLP, ...
    kL, clp, cle, f1, f2, h1, h2, w, l, shat, deltaS, eta, k2, ...
    hs, tg1, ee0, p, deltaF, V, deltaL, deltaE, deltaP, ...
    n] = Franken_params ;
    
% rename variables
m0 = Y(1);
m1 = Y(2);
m2 = Y(3);
k0 = Y(4);
ka = Y(5);
e0 = Y(6);
ef = Y(7);
ea = Y(8);
Ra = Y(9);
Rb = Y(10);

FSHp = Y(11);
FSH = Y(12);
LHp = Y(13);
LH = Y(14);
PH = Y(15);
OM = Y(16);
LB = Y(17);
S = Y(18);
E2 = Y(19);
P4 = Y(20);

%% ode system %%

% dimers using QSS (scaled)
Da = calphaP / calphaM .* Ra.^2 * K2; % homodimer
Dab = cbetaP / cbetaM .* K2 .* Ra .* Rb ; % heterodimer

% macrophage
dm0dt = mum + etaM .* (1- (m0 + m1 + m2)./mc).*(m0 + m1 + m2) - thetaM .* (ka / (cka + ka)).* m0 - beta1 .* (ef + ea) .* m0 - beta2 .* m0 - deltaM0 .* m0;
dm1dt = beta1 .* (ef + ea).* m0 + thetaM.*(ka ./ (cka + ka)).* m0 + beta21.* (1 + (phi4 .* Dab.^n4 ./ (Dab.^n4 + H4.^n4))).* m2 - beta12.*(1 + (phi3 .* Da.^n3 ./ (Da.^n3 + H3.^n3))) .* m1 - deltaM1 .* m1;
dm2dt = beta2 .* m0 + beta12 .* (1 + (phi3 .* Da.^n3 ./ (Da.^n3 + H3.^n3))).* m1 - beta21 .*(1 + (phi4 .* Dab.^n4 ./ (Dab.^n4 + H4.^n4))) .* m2 - deltaM2 .* m2;

% nk cells
dk0dt = (1 + (phi1 .* E2)) .* muk - thetaK .* k0 .* (m1 ./ (cm1 + m1)) - deltaK0 .* k0;
dkadt = thetaK .* k0 .* (m1 ./ (cm1 + m1)) + etaK .* (1 - (ka + k0)./kc).* ka - sigma .* (ef + ea).* ka - deltaKa * ka;

% endometrial cells
de0dt = mue .* P4 - deltaE0 .* (rho0 .* e0 + (1-rho0).* e0);
defdt = deltaE0 .* rho0 .* e0 - (phi2 .* H2.^n2 ./ (Rb.^n2 + H2.^n2)) .* omega .*(gamma .* ka + (1-gamma).* m1).* ef - rhoF .* ef - deltaEf .* ef;
deadt = rhoF .* ef + etaE .* ea .* (1 - (ea ./ ec)) .* (m2 ./ (cm2 + m2)) - (phi5 .* H5.^n5 ./ (Rb.^n5 + H5.^n5)) .* omega.* (gamma.* ka + (1-gamma).* m1).*ea - deltaEa .* ea;

% receptor dynamics
%dRadt = kalphaP .* E2 .* (K2 - Dab - Da - Ra) - kalphaM .* Ra - cbetaP .* Ra .* Rb + cbetaM .* Dab - calphaP .* Ra.^2 + calphaM .* Da; % bound ERa
%dRbdt = kbetaP .* E2 .* (K1 - Dab - Rb) - kbetaM .* Rb - cbetaP .* Ra .* Rb + cbetaM .* Dab; % bound ERb

% scaled receptor dynamics
dRadt = kalphaP .* E2 .* (1 - Da - Ra - (K1/K2) .* Dab) - kalphaM .* Ra - cbetaP .* Ra .* Rb .* K1 + cbetaM .* K1/K2 .* Dab  - calphaP .* Ra.^ 2 .* K2 + calphaM .* Da ;
dRbdt = kbetaP .* E2 .* (1 - Dab - Rb) - kbetaM .* Rb - cbetaP .* K2 .* Ra .* Rb + cbetaM .* Dab ;

% FSH and LH
dFSHpdt = (vf ./ (1 + cfi .* ((S .* LB)./ (KiFI + S*LB)))) - kF .* ((1 + cfp .* P4) ./(1 + cfe .* E2.^2)).* FSHp ;
dFSHdt = 1/V .* kF .* ((1 + cfp .* P4) ./ (1 + cfe .* E2 .^2)).* FSHp - deltaF .* FSH ;
dLHpdt = v0L + ((v1L .* E2 .^n) ./ (KmL .^n + E2 .^n)) .* 1 / (1 + P4 ./ KiLP) - kL .* ((1 + clp .* P4)./ (1 + cle .* E2)) .* LHp ;
dLHdt = 1/V .* kL .* ((1 + clp .* P4) ./ (1 + cle .* E2)) .* LHp - deltaL .*LH ;

% Follicule Dynamics
dPHdt = PH .* ((f1 .* FSH .^2) ./ (h1 .^2 + FSH .^ 2) - ((f2 .* LH .^ 2)./ (h2 .^2 + LH .^ 2)));
dOMdt = ((f2 .* LH .^2) ./ (h2 .^2 + LH .^2)) .* PH - w .* S .* OM;
dLBdt = w .* S .* OM - l .* (1-S) .* LB;
dSdt = shat .* (LH .^4) ./ (LH .^ 4 + hs .^ 4) .* (1 - S) - deltaS .* S;

% Estrogen and Progesterone
dE2dt = ee0 - deltaE .* E2 + tg1 .* (LH ./ (LH + k2)) .* (PH + eta .* LB .* S) + psi .* ea + kbetaM .* Rb - kbetaP .* E2 .* (1 - Dab - Rb) + kalphaM .* K2/K1 .* Ra - kalphaP .* E2 .* (K2/K1 - Dab - K2/K1 .* Da - K2/K1 .* Ra);
dP4dt = -deltaP .* P4 + p .* LB .* S;

%% output %%

out = [dm0dt; dm1dt; dm2dt; dk0dt; dkadt; de0dt; defdt; deadt; ...
    dRadt; dRbdt; dFSHpdt ; dFSHdt ; dLHpdt ;dLHdt ; dPHdt ; dOMdt ; dLBdt ; dSdt ; dE2dt ; dP4dt];


