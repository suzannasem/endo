%%%% Reproducing Miller(2025) simulations %%%%

%%% function file %%%

function out = Miller_fn(t,Y,params, beta1, omega, rho0)

%% Function based on the work by Claire Miller et. al (2025)
% inputs: params (structure) - holds all parameter values for system, minus
% beta1, omega, and rho0: input parameters dictating immune function and
% retrograde influx rate

%% parameters %%

% macrophage params
mum = params.mum ;
etaM = params.etaM ;
mc = params.mc ;
thetaM = params.thetaM ;
cka = params.cka ;
beta2 = params.beta2 ;
deltaM0 = params.deltaM0 ;
beta21 = params.beta21 ;
beta12 = params.beta12 ;
deltaM1 = params.deltaM1 ;
deltaM2 = params.deltaM2 ;

% nk params
muk = params.muk ;
thetaK = params.thetaK ;
cm1 = params.cm1 ;
deltaK0 = params.deltaK0 ;
deltaKa = params.deltaK0 ;
etaK = params.etaK ;
kc = params.kc ;
sigma = params.sigma ;

% endo params
ae = 2e4;
be=160;
de=12;

mue = ae.*(sin(pi.*(t+de)./28)).^be;  
%mue = params.mue ; % only use in surrogate model
deltaE0 = params.deltaE0 ; 
gamma = params.gamma ;
rhoF = params.rhoF ;
deltaEf = params.deltaEf ;
etaE = params.etaE ;
ec = params.ec ;
cm2 = params.cm2 ;
deltaEa = params.deltaEa ;

% rename variable
m0 = Y(1);
m1 = Y(2);
m2 = Y(3);
k0 = Y(4);
ka = Y(5);
e0 = Y(6);
ef = Y(7);
ea = Y(8);

%% ode system %%

% macrophage
dm0dt = mum + etaM .* (1- (m0 + m1 + m2)./mc).*(m0 + m1 + m2) - thetaM .* (ka / (cka + ka)).* m0 - beta1 .* (ef + ea) .* m0 - beta2 .* m0 - deltaM0 .* m0;
dm1dt = beta1 .* (ef + ea).* m0 + thetaM.*(ka ./ (cka + ka)).* m0 + beta21.* m2 - beta12.* m1 - deltaM1 .* m1;
dm2dt = beta2 .* m0 + beta12 .* m1 - beta21 .* m2 - deltaM2 .* m2;

% nk cells
dk0dt = muk - thetaK .* k0 .* (m1 ./ (cm1 + m1)) - deltaK0 .* k0;
dkadt = thetaK .* k0 .* (m1 ./ (cm1 + m1)) + etaK .* (1 - (ka + k0)./kc).* ka - sigma .* (ef + ea).* ka - deltaKa * ka;

% endometrial cells
de0dt = mue - deltaE0 .* (rho0 .* e0 + (1-rho0).* e0);
defdt = deltaE0 .* rho0 .* e0 - omega .*(gamma .* ka + (1-gamma).* m1).* ef - rhoF .* ef - deltaEf .* ef;
deadt = rhoF .* ef + etaE .* ea .* (1 - (ea ./ ec)) .* (m2 ./ (cm2 + m2)) - omega.* (gamma.* ka + (1-gamma).* m1).*ea - deltaEa .* ea;

%% output %%
out = [dm0dt; dm1dt; dm2dt; dk0dt; dkadt; de0dt; defdt; deadt];
