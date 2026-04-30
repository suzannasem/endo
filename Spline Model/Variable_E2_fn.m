%%%% Incorporating E2 and P4 %%%%

%%% function file %%%

function out = Variable_E2_fn(t,Y, beta1, omega, E2, P4, K1, K2)
%% parameters %%

[mum,muk,etaE,etaK,etaM,deltaM0,deltaM1,deltaM2,deltaK0,deltaKa,deltaE0,deltaEa, deltaEf,sigma,mc, ...
              kc,ec,cka,cm1,cm2,beta2,beta12,beta21,thetaM,thetaK,gamma, ...
                 rhoF, rho0, mue, ...
                 phi1, n1, H1, phi2, n2, H2, ...
                 kalphaM, kalphaP, kbetaM, kbetaP, calphaM, calphaP, cbetaM, cbetaP, ... % K1, K2, ...
                 phi3, n3, H3, phi4, n4, H4, ...
                 phi5, n5, H5] = Constant_E2_params;
    
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

%% ode system %%

% E2 and P4
E2_now = E2(t);
P4_now = P4(t);

% cyclic (sine curve) as influx
%ae = 2e4;
%be=160;
%de=12;
% mue = ae.*(sin(pi.*(t+de)./28)).^be;  
% mue = Constant_E2_params.mue ; % turn to make heatmaps run faster

% dimers using QSS
%Da = calphaP / calphaM .* Ra.^2 ; % homodimer
%Dab = cbetaP / cbetaM .* Ra .* Rb ; % heterodimer
% scaled version
Da = calphaP / calphaM .* Ra.^2 * K2; % homodimer
Dab = cbetaP / cbetaM .* K2 .* Ra .* Rb ; % heterodimer

% macrophage
dm0dt = mum + etaM .* (1- (m0 + m1 + m2)./mc).*(m0 + m1 + m2) - thetaM .* (ka / (cka + ka)).* m0 - beta1 .* (ef + ea) .* m0 - beta2 .* m0 - deltaM0 .* m0;
dm1dt = beta1 .* (ef + ea).* m0 + thetaM.*(ka ./ (cka + ka)).* m0 + beta21.* (1 + (phi4 .* Dab.^n4 ./ (Dab.^n4 + H4.^n4))).* m2 - beta12.*(1 + (phi3 .* Da.^n3 ./ (Da.^n3 + H3.^n3))) .* m1 - deltaM1 .* m1;
dm2dt = beta2 .* m0 + beta12 .* (1 + (phi3 .* Da.^n3 ./ (Da.^n3 + H3.^n3))).* m1 - beta21 .*(1 + (phi4 .* Dab.^n4 ./ (Dab.^n4 + H4.^n4))) .* m2 - deltaM2 .* m2;

% nk cells
dk0dt = (1 + (phi1 .* E2_now.^n1 ./ (E2_now.^n1 + H1.^n1))).*muk - thetaK .* k0 .* (m1 ./ (cm1 + m1)) - deltaK0 .* k0;
dkadt = thetaK .* k0 .* (m1 ./ (cm1 + m1)) + etaK .* (1 - (ka + k0)./kc).* ka - sigma .* (ef + ea).* ka - deltaKa * ka;

% endometrial cells
de0dt = mue .* P4_now - deltaE0 .* (rho0 .* e0 + (1-rho0).* e0);
defdt = deltaE0 .* rho0 .* e0 - (phi2 .* H2.^n2 ./ (Rb.^n2 + H2.^n2)) .* omega .*(gamma .* ka + (1-gamma).* m1).* ef - rhoF .* ef - deltaEf .* ef;
deadt = rhoF .* ef + etaE .* ea .* (1 - (ea ./ ec)) .* (m2 ./ (cm2 + m2)) - (phi5 .* H5.^n5 ./ (Rb.^n5 + H5.^n5)) .* omega.* (gamma.* ka + (1-gamma).* m1).*ea - deltaEa .* ea;

% receptor dynamics
%dRadt = kalphaP .* E2_now .* (K2 - Dab - Da - Ra) - kalphaM .* Ra - cbetaP .* Ra .* Rb + cbetaM .* Dab - calphaP .* Ra.^2 + calphaM .* Da; % bound ERa
%dRbdt = kbetaP .* E2_now .* (K1 - Dab - Rb) - kbetaM .* Rb - cbetaP .* Ra .* Rb + cbetaM .* Dab; % bound ERb
% scaled version
dRadt = kalphaP .* E2_now .* (1 - Da - Ra - (K1/K2) .* Dab) - kalphaM .* Ra - cbetaP .* Ra .* Rb .* K1 + cbetaM .* K1/K2 .* Dab  - calphaP .* Ra.^ 2 .* K2 + calphaM .* Da ;
dRbdt = kbetaP .* E2_now .* (1 - Dab - Rb) - kbetaM .* Rb - cbetaP .* K2 .* Ra .* Rb + cbetaM .* Dab ;


%% output %%

out = [dm0dt; dm1dt; dm2dt; dk0dt; dkadt; de0dt; defdt; deadt; ...
    dRadt; dRbdt];
