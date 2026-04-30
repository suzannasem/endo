%%%% LESIOn: Graham (2023) + Miller (2025) + Semaan (2026) + local estradiol %%%%

function out = LESIOn_fn(t,Y, K1, K2, psi, p)
   
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
LE2 = Y(21);

%% ode system %%

% dimers using QSS
Da = p.calphaP / p.calphaM .* Ra.^2 * K2; % homodimer
Dab = p.cbetaP / p.cbetaM .* K2 .* Ra .* Rb ; % heterodimer

% macrophage
dm0dt = p.mum + p.etaM .* (1- (m0 + m1 + m2)./p.mc).*(m0 + m1 + m2) - p.thetaM .* (ka / (p.cka + ka)).* m0 - p.beta1 .* (ef + ea) .* m0 - p.beta2 .* m0 - p.deltaM0 .* m0;
dm1dt = p.beta1 .* (ef + ea).* m0 + p.thetaM.*(ka ./ (p.cka + ka)).* m0 + p.beta21.* (1 + (p.phi4 .* Dab.^p.n4 ./ (Dab.^p.n4 + p.H4.^p.n4))).* m2 - p.beta12.*(1 + (p.phi3 .* Da.^p.n3 ./ (Da.^p.n3 + p.H3.^p.n3))) .* m1 - p.deltaM1 .* m1;
dm2dt = p.beta2 .* m0 + p.beta12 .* (1 + (p.phi3 .* Da.^p.n3 ./ (Da.^p.n3 + p.H3.^p.n3))).* m1 - p.beta21 .*(1 + (p.phi4 .* Dab.^p.n4 ./ (Dab.^p.n4 + p.H4.^p.n4))) .* m2 - p.deltaM2 .* m2;

% nk cells
dk0dt = (1 + (p.phi1 .* E2)) .* p.muk - p.thetaK .* k0 .* (m1 ./ (p.cm1 + m1)) - p.deltaK0 .* k0;
dkadt = p.thetaK .* k0 .* (m1 ./ (p.cm1 + m1)) + p.etaK .* (1 - (ka + k0)./p.kc).* ka - p.sigma .* (ef + ea).* ka - p.deltaKa * ka;

% endometrial cells
de0dt = p.mue .* P4 - p.deltaE0 .* (p.rho0 .* e0 + (1-p.rho0).* e0);
defdt = p.deltaE0 .* p.rho0 .* e0 - (p.phi2 .* p.H2.^p.n2 ./ (Rb.^p.n2 + p.H2.^p.n2)) .* p.omega .*(p.gamma .* ka + (1-p.gamma).* m1).* ef - p.rhoF .* ef - p.deltaEf .* ef;
deadt = p.rhoF .* ef + p.etaE .* ea .* (1 - (ea ./ p.ec)) .* (m2 ./ (p.cm2 + m2)) - (p.phi5 .* p.H5.^p.n5 ./ (Rb.^p.n5 + p.H5.^p.n5)) .* p.omega.* (p.gamma.* ka + (1-p.gamma).* m1).*ea - p.deltaEa .* ea;

% receptor dynamics
%dRadt = kalphaP .* E2 .* (K2 - Dab - Da - Ra) - kalphaM .* Ra - cbetaP .* Ra .* Rb + cbetaM .* Dab - calphaP .* Ra.^2 + calphaM .* Da; % bound ERa
%dRbdt = kbetaP .* E2 .* (K1 - Dab - Rb) - kbetaM .* Rb - cbetaP .* Ra .* Rb + cbetaM .* Dab; % bound ERb

% scaled receptor dynamics
dRadt = p.kalphaP .* E2 .* (1 - Da - Ra - (K1/K2) .* Dab) - p.kalphaM .* Ra - p.cbetaP .* Ra .* Rb .* K1 + p.cbetaM .* K1/K2 .* Dab  - p.calphaP .* Ra.^ 2 .* K2 + p.calphaM .* Da ;
dRbdt = p.kbetaP .* E2 .* (1 - Dab - Rb) - p.kbetaM .* Rb - p.cbetaP .* K2 .* Ra .* Rb + p.cbetaM .* Dab ;

% FSH and LH
dFSHpdt = (p.vf ./ (1 + p.cfi .* ((S .* LB)./ (p.KiFI + S*LB)))) - p.kF .* ((1 + p.cfp .* P4) ./(1 + p.cfe .* E2.^2)).* FSHp ;
dFSHdt = 1/p.V .* p.kF .* ((1 + p.cfp .* P4) ./ (1 + p.cfe .* E2 .^2)).* FSHp - p.deltaF .* FSH ;
dLHpdt = p.v0L + ((p.v1L .* E2 .^p.n) ./ (p.KmL .^p.n + E2 .^p.n)) .* 1 / (1 + P4 ./ p.KiLP) - p.kL .* ((1 + p.clp .* P4)./ (1 + p.cle .* E2)) .* LHp ;
dLHdt = 1/p.V .* p.kL .* ((1 + p.clp .* P4) ./ (1 + p.cle .* E2)) .* LHp - p.deltaL .*LH ;

% Follicule Dynamics
dPHdt = PH .* ((p.f1 .* FSH .^2) ./ (p.h1 .^2 + FSH .^ 2) .* (1 ./ (p.pk .* LE2 + 1)) - ((p.f2 .* LH .^ 2)./ (p.h2 .^2 + LH .^ 2)));
dOMdt = ((p.f2 .* LH .^2) ./ (p.h2 .^2 + LH .^2)) .* PH - p.w .* S .* OM;
dLBdt = p.w .* S .* OM - p.l .* (1-S) .* LB;
dSdt = p.shat .* (LH .^4) ./ (LH .^ 4 + p.hs .^ 4) .* (1 - S) - p.deltaS .* S;

% Estrogen and Progesterone
dE2dt = p.ee0 - p.deltaE .* E2 + p.tg1 .* (LH ./ (LH + p.k2)) .* (PH + p.eta .* LB .* S) + p.kbetaM .* Rb - p.kbetaP .* E2 .* (1 - Dab - Rb) + p.kalphaM .* K2/K1 .* Ra - p.kalphaP .* E2 .* (K2/K1 - Dab - K2/K1 .* Da - K2/K1 .* Ra);
dP4dt = -p.deltaP .* P4 + p.p .* LB .* S;

% Local Estrogen
dLE2dt = psi .* ea - p.deltaLE2 .* LE2 ;

%% output %%

out = [dm0dt; dm1dt; dm2dt; dk0dt; dkadt; de0dt; defdt; deadt; ...
    dRadt; dRbdt; dFSHpdt ; dFSHdt ; dLHpdt ;dLHdt ; dPHdt ; dOMdt ; dLBdt ; dSdt ; dE2dt ; dP4dt ; ...
    dLE2dt]; 


