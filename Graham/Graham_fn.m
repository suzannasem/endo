%% Function File for Graham Reduced Model %%

function out = Graham(t, Y, params)

%% parameters %%

[vf, KiFI, kF, cfp, cfe, cfi, v0L, v1L, KmL, KiLP, ...
    kL, clp, cle, f1, f2, h1, h2, w, l, shat, deltaS, eta, k2, ...
    hs, tg1, e0, p, deltaF, V, deltaL, deltaE, deltaP, ...
    n] = graham_params ;

% rename variable
FSHp = Y(1);
FSH = Y(2);
LHp = Y(3);
LH = Y(4);
PH = Y(5);
OM = Y(6);
LB = Y(7);
S = Y(8);
E2 = Y(9);
P4 = Y(10);

%% ode system %%

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
dE2dt = e0 - deltaE .* E2 + tg1 .* (LH ./ (LH + k2)) .* (PH + eta .* LB .* S);
dP4dt = -deltaP .* P4 + p .* LB .* S;

%% output %%
out = [dFSHpdt ; dFSHdt ; dLHpdt ;dLHdt ; dPHdt ; dOMdt ; dLBdt ; dSdt ; dE2dt ; dP4dt];

