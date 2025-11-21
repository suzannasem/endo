%%%% Reproducing Miller(2025) simulations %%%%
%%% Simulation file %%%

%% Define parameters %%

% macrophage params
params.mum = 10^4; % cells / mL * day
params.etaM = 0.2; % 1/day
params.mc = 10^6; % cells/mL
params.thetaM = 5.0 * 10^(-8); % 1/day
params.cka = 1e5 ; %10^6; % cells/mL %1e5 new
params.beta2 = 10^(-3); % 1/day
params.deltaM0 = 0.02; %1/day
params.beta21 = 5.0 * 10^(-5); % 1/day
params.beta12 = 5.0 * 10^(-5); % 1/day
params.deltaM1 = 0.02 ; % 1/day
params.deltaM2 = 0.02; % 1/day

% nk params
params.muk = 3.5e3 ; % 10^4; % cells / mL * day
params.thetaK = 0.7 ; %0.2; % 1/day
params.cm1 = 10^5; % cells/mL
params.deltaK0 = 0.02; % 1/day
params.etaK = 0.02; % 1/day
params.kc = 10^6 ; % cells / mL
params.sigma = 10^(-5); % mL / cells * day

% endo params
%params.mue = 1260; % cells / mL * day - avg
params.deltaE0 = 1.0 ; % 1/day
params.gamma = 0.8 ;
params.rhoF = 0.1; % 1/day
params.deltaEf = 0.14 ; % 1/day
params.etaE = 1.0; % 1/day
params.ec = 10^9 ; % cells/mL
params.cm2 = 100; % cells/mL
params.deltaEa = 0.14; % 1/day 


%% Initial conditions %%

m00 = 9.0 * 10^5 ;
m10 = 100 ;
m20 = 4.5 * 10^4 ;
k00 = 5 * 10^5 ; 
ka0 = 2 * 10^4 ; 
e00 = 0 ; ef0 = 0; ea0 = 0;

init = [m00, m10, m20, k00, ka0, e00, ef0, ea0] ;
% alternate initial conditions
%y0 = [898042.5; 9637.075; 44814.18; 42929.06; 296376.2; 1260; 105.3528; 107.2753];
%dY = [-1; 1; -1; -1; -1; 0; 1; 1];
%y0 = y0.*(1 + 0.1*dY);

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

%% ODE Simulations %%

tspan = 0:0.5:1200 ;
[t,Y] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6), 10^(-5), 0.1),tspan,init, opts); % given params - 1
[t1,Y1] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6.9), 10^(-4.9), 0.1),tspan,init, opts); % scenario 2 on fig 4 b
[t2,Y2] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6.4), 10^(-5.3), 0.1),tspan,init, opts); % 3
[t11,Y11] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6.25), 10^(-5.28), 0.1),tspan,init, opts); % 4
[t22,Y22] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6.0), 10^(-5.7), 0.1),tspan, init, opts); % 5

% clearance rates (rhoF)
[t3,Y3] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6), 10^(-5), 0),tspan,init); % no retrograde menstruation
[t4,Y4] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6), 10^(-5), 0.005),tspan,init); % 0.005
[t5,Y5] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6), 10^(-5), 0.01),tspan,init); % 0.01
[t6,Y6] = ode15s(@(t,y) Miller_fn(t,y,params, 10^(-6), 10^(-5), 0.2),tspan,init); % lots of influx

% default (scenario 1)
m0 = Y(:,1);
m1 = Y(:,2);
m2 = Y(:,3);
k0 = Y(:,4);
ka = Y(:,5);
e0 = Y(:,6);
ef = Y(:,7);
ea = Y(:,8);

% scenario 2
m01 = Y1(:,1);
m11 = Y1(:,2);
m21 = Y1(:,3);
k01 = Y1(:,4);
ka1 = Y1(:,5);
e01 = Y1(:,6);
ef1 = Y1(:,7);
ea1 = Y1(:,8);

% scenario 3
m02 = Y2(:,1);
m12 = Y2(:,2);
m22 = Y2(:,3);
k02 = Y2(:,4);
ka2 = Y2(:,5);
e02 = Y2(:,6);
ef2 = Y2(:,7);
ea2 = Y2(:,8);

% scenario 4
m011 = Y11(:,1);
m111 = Y11(:,2);
m211 = Y11(:,3);
k011 = Y11(:,4);
ka11 = Y11(:,5);
e011 = Y11(:,6);
ef11 = Y11(:,7);
ea11 = Y11(:,8);

% scenario 5
m022 = Y22(:,1);
m122 = Y22(:,2);
m222 = Y22(:,3);
k022 = Y22(:,4);
ka22 = Y22(:,5);
e022 = Y22(:,6);
ef22 = Y22(:,7);
ea22 = Y22(:,8);

% clearance rates

% rhoF = 0
m03 = Y3(:,1);
m13 = Y3(:,2);
m23 = Y3(:,3);
k03 = Y3(:,4);
ka3 = Y3(:,5);
e03 = Y3(:,6);
ef3 = Y3(:,7);
ea3 = Y3(:,8);

% rhoF = 0.005
m04 = Y4(:,1);
m14 = Y4(:,2);
m24 = Y4(:,3);
k04 = Y4(:,4);
ka4 = Y4(:,5);
e04 = Y4(:,6);
ef4 = Y4(:,7);
ea4 = Y4(:,8);

% rhoF = 0.01
m05 = Y5(:,1);
m15 = Y5(:,2);
m25 = Y5(:,3);
k05 = Y5(:,4);
ka5 = Y5(:,5);
e05 = Y5(:,6);
ef5 = Y5(:,7);
ea5 = Y5(:,8);

% rhoF = 0.2
m06 = Y6(:,1);
m16 = Y6(:,2);
m26 = Y6(:,3);
k06 = Y6(:,4);
ka6 = Y6(:,5);
e06 = Y6(:,6);
ef6 = Y6(:,7);
ea6 = Y6(:,8);

%% Plotting %%

% Figure 1

figure(1)

subplot(2,2,1)
plot(t, k0, '--', t, ka, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1000,1080])
ylim([0,400e3])
title("Natural Killer Cells")

subplot(2,2,2)
plot(t, e0, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1000,1080])
ylim([0, 15e3])
title("Eutopic Endometrial Cells")

subplot(2,2,3)
semilogy(t, m0, '--', t, m1, '-', t, m2, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1000,1080])
title("Macrophage")

subplot(2,2,4)
plot(t, ea, '-', t, ef, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1000,1080])
ylim([0,0.7e3])
title("Ectopic Endometrial Cells")

% different disease states
figure(2)

subplot(3,1,1)
semilogy(t2, ea2, '-.', t11, ea11, '.', t22, ea22, ':', 'LineWidth', 2)
legend('3','4','5')
xlim([1000,1080])
ylabel("E_A (cells/mL)")
ylim([10^0, 10^9.1])
title("Attached Cells in Disease State")

subplot(3,1,2)
plot(t, ea, '-', t1, ea1, '--', 'LineWidth', 2)
legend('1','2')
xlim([1000,1080])
ylabel("E_A (cells/mL)")
title("Attached Cells in Low/No Disease State")

subplot(3,1,3)
plot(t, ef, '-', t1, ef1, '--', t2, ef2, '-.', t11, ef11, '.', t22, ef22, ':', 'LineWidth', 2)
legend('1','2','3','4','5')
xlim([1000,1080])
ylabel("E_F (cells/mL)")
title("Endometrial Cells in Peritoneal Fluid")


% varying retrograde influx
figure(3)

subplot(3,2,1)
axis off
hold on;
h1 = plot(nan, nan, 'b--', 'DisplayName', '0');
h2 = plot(nan, nan, 'r-.', 'DisplayName', '0.005');
h3 = plot(nan, nan, 'y.', 'DisplayName', '0.01');
h4 = plot(nan, nan, 'v--', 'DisplayName', '0.1');
h5 = plot(nan, nan, 'g:', 'DisplayName', '0.2');
hold off
legend([h1 h2 h3 h4 h5], 'Location', 'best', 'Box', 'on');
title('\rho Legend')

subplot(3,2,2)
plot(t3, ka3, '--', t4, ka4, '-.', t5, ka5, '.', t, ka, '-', t6, ka6, ':', 'LineWidth', 2)
%legend("0", "0.005", "0.01", "0.1", "0.2")
xlabel("Time in days")
ylabel("K_A Cells / mL (\times 10^3)")
xlim([1000, 1080])
ylim([0,300e3])
title("Activated NK cells, Varying Retrograde Clearance \rho_0")

subplot(3,2,3)
plot(t3, ef3, '--', t4, ef4, '-.', t5, ef5, '.', t, ef, '-', t6, ef6, ':', 'LineWidth', 2)
%legend("0", "0.005", "0.01", "0.1", "0.2")
xlabel("Time in days")
ylabel("E_f Cells / mL")
xlim([1000, 1080])
%ylim([0,1e3])
title("Endometrial Cells in Peritoneal Fluid, Varying Retrograde Clearance \rho_0")

m1_p = m1./(m0 + m1 + m2) ;
m13_p = m13./(m03 + m13 + m23) ;
m14_p = m14./(m04 + m14 + m24) ;
m15_p = m15./(m04 + m15 + m25) ;
m16_p = m16./(m06 + m16 + m26) ;

subplot(3,2,4)
plot(t3, m13_p, '--', t4, m14_p, '-.', t5, m15_p, '.', t, m1_p, '-', t6, m16_p, ':', 'LineWidth', 2)
%legend("0", "0.005", "0.01", "0.1", "0.2")
xlabel("Time in days")
ylabel("Prop of M1 cells")
xlim([1000, 1080])
title("Proportion of Inflammatory-type Macrophage")

m2_p = m2./(m0 + m1 + m2) ;
m23_p = m23./(m03 + m13 + m23) ;
m24_p = m24./(m04 + m14 + m24) ;
m25_p = m25./(m04 + m15 + m25) ;
m26_p = m26./(m06 + m16 + m26) ;

subplot(3,2,5)
plot(t3, m23_p, '--', t4, m24_p, '-.', t5, m25_p, '.', t, m2_p, '-', t6, m26_p, ':', 'LineWidth', 2)
%legend("0", "0.005", "0.01", "0.1", "0.2")
xlabel("Time in days")
ylabel("Prop of M2 cells")
xlim([1000,1080])
ylim([0.046, 0.049])
title("Proportion of Non-Inflammatory-type Macrophage")

subplot(3,2,6)
plot(t3, ea3, '--', t4, ea4, '-.', t5, ea5, '.', t, ea, '-', t6, ea6, ':', 'LineWidth', 2)
%legend("0", "0.005", "0.01", "0.1", "0.2")
xlabel("Time in days")
ylabel("E_f Cells / mL")
xlim([1000, 1080])
%ylim([0,1e3])
title("Lesion Cells")

%% heatmaps from fig 4 in text

% sequences of parameter values to test
omega_list = linspace(1e-6,1e-4,5);
beta1_list = linspace(1e-7,1e-5,5);

% empty matrix for storing ea
ea_vals_low = zeros(5,5) ;
ea_vals_mod = zeros(5,5) ;
ea_vals_high = zeros(5,5) ;

% empty matrix for storing m0
m0_vals_low = zeros(5,5) ;
m0_vals_mod = zeros(5,5) ;
m0_vals_high = zeros(5,5) ;

% empty matrix for storing m1
m1_vals_low = zeros(5,5) ;
m1_vals_mod = zeros(5,5) ;
m1_vals_high = zeros(5,5) ;

% empty matrix for storing m2
m2_vals_low = zeros(5,5) ;
m2_vals_mod = zeros(5,5) ;
m2_vals_high = zeros(5,5) ;

% timing the simulation
tic

% loop to fill matrix
for j = 1:length(beta1_list)
    for i = 1:length(omega_list)
        % run simulations
        [t_ea_low,Y_ea_low] = ode15s(@(t,y) Miller_fn(t,y,params, beta1_list(j), omega_list(i), 0.01),tspan,init, opts);
        [t_ea_mod,Y_ea_mod] = ode15s(@(t,y) Miller_fn(t,y,params, beta1_list(j), omega_list(i), 0.1),tspan,init, opts);
        [t_ea_high,Y_ea_high] = ode15s(@(t,y) Miller_fn(t,y,params, beta1_list(j), omega_list(i), 0.2),tspan,init, opts);
        % storing ea values for heatmap 1
        ea_vals_low(j,i) = Y_ea_low(1020,8);
        ea_vals_mod(j,i) = Y_ea_mod(1020,8);
        ea_vals_high(j,i) = Y_ea_high(1020,8);
        % storing m0 values for proportions
        m0_vals_low(j,i) = Y_ea_low(1020,1);
        m0_vals_mod(j,i) = Y_ea_mod(1020,1);
        m0_vals_high(j,i) = Y_ea_high(1020,1);
        % storing m1 values for heatmaps 2
        m1_vals_low(j,i) = Y_ea_low(1020,2);
        m1_vals_mod(j,i) = Y_ea_mod(1020,2);
        m1_vals_high(j,i) = Y_ea_high(1020,2);
        % storing m2 values for heatmaps 3
        m2_vals_low(j,i) = Y_ea_low(1020,3);
        m2_vals_mod(j,i) = Y_ea_mod(1020,3);
        m2_vals_high(j,i) = Y_ea_high(1020,3);
    end
end

% timing
elapsedTime = toc;  % stop timer
fprintf('Loop runtime: %.2f seconds\n', elapsedTime);

% plot heatmap - Figure 4(b) in text
figure(4)

subplot(1,3,1)
imagesc([1e-6,1e-4],[1e-7,1e-5],log10(ea_vals_low))
set(gca, 'XScale', 'log', 'YScale', 'log','YDir','normal')
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('Low influx: \rho_0 = 0.01')
cb = colorbar ;
cb.Ticks = -3:9 ;
cb.Label.String = 'Lesion Cells, in Powers of 10' ;

subplot(1,3,2)
imagesc([1e-6,1e-4],[1e-7,1e-5],log10(ea_vals_mod))
set(gca, 'XScale', 'log', 'YScale', 'log','YDir','normal')
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('Moderate influx: \rho_0 = 0.1')
cb = colorbar ;
cb.Ticks = -3:9 ;
cb.Label.String = 'Lesion Cells, in Powers of 10' ;

subplot(1,3,3)
imagesc([1e-6,1e-4],[1e-7,1e-5],log10(ea_vals_high))
set(gca, 'XScale', 'log', 'YScale', 'log','YDir','normal')
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('High influx: \rho_0 = 0.2')
cb = colorbar ;
cb.Ticks = -3:9 ;
cb.Label.String = 'Lesion Cells, in Powers of 10' ;

% proportions of m1 and m2 macrophage
m1_prop_low = m1_vals_low ./ (m0_vals_low + m1_vals_low + m2_vals_low) ;
m2_prop_low = m2_vals_low ./ (m0_vals_low + m1_vals_low + m2_vals_low) ;

m1_prop_mod = m1_vals_mod ./ (m0_vals_mod + m1_vals_mod + m2_vals_mod) ;
m2_prop_mod = m2_vals_mod ./ (m0_vals_mod + m1_vals_mod + m2_vals_mod) ;

m1_prop_high = m1_vals_high ./ (m0_vals_high + m1_vals_high + m2_vals_high) ;
m2_prop_high = m2_vals_high ./ (m0_vals_high + m1_vals_high + m2_vals_high) ;

% m1 & m2 proportions
m11_p = m11./(m01 + m11 + m21) ; % scenario 2
m12_p = m12./(m02 + m12 + m22) ; % scenario 3
m111_p = m111./(m011 + m111 + m211) ; % scenario 4
m122_p = m122./(m022 + m122 + m222) ; % scenario 5

m21_p = m21./(m01 + m11 + m21) ; % scenario 2
m22_p = m22./(m02 + m12 + m22) ; % scenario 3
m211_p = m211./(m011 + m111 + m211) ; % scenario 4
m222_p = m222./(m022 + m122 + m222) ; % scenario 5

% plot heatmap - Figure 5a/b in text
figure(5)

subplot(2,3,1)
heatmap(omega_list, beta1_list, m1_prop_low) ;
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('Low influx: \rho_0 = 0.01')
colorbar 

subplot(2,3,2)
heatmap(omega_list, beta1_list, m1_prop_mod)
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('Moderate influx: \rho_0 = 0.1')
colorbar

subplot(2,3,3)
heatmap(omega_list, beta1_list, m1_prop_high)
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('High influx: \rho_0 = 0.2')
colorbar
cb.Ticks = -3:9 ;

subplot(2,3,4)
heatmap(omega_list, beta1_list, m2_prop_low)
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('Low influx: \rho_0 = 0.01')
colorbar

subplot(2,3,5)
heatmap(omega_list, beta1_list, m2_prop_mod)
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('Moderate influx: \rho_0 = 0.1')
clim([1e-4, 1e-1]);
colorbar

subplot(2,3,6)
heatmap(omega_list, beta1_list, m2_prop_high)
xlabel('Clearance rate, \omega (c^{-1}\cdot mL \cdot d^{-1})')
ylabel('Detection Rate, \beta_1 (c^{-1}\cdot mL \cdot d^{-1})')
title('High influx: \rho_0 = 0.2')
colorbar

figure(6)

subplot(1,2,1)
semilogy(t, m1_p,'--', t1, m11_p, '-', t2, m12_p, ':', t11, m111_p, '-.', t22, m122_p, '.', 'LineWidth', 2)
legend("1", "2", "3", "4", "5")
xlabel("Time in days")
ylabel("M_1 Proportion")
xlim([1010, 1090])

subplot(1,2,2)
semilogy(t, m2_p,'--', t1, m21_p, '-', t2, m22_p, ':', t11, m211_p, '-.', t22, m222_p, '.', 'LineWidth', 2)
legend("1", "2", "3", "4", "5")
xlabel("Time in days")
ylabel("M_2 Proportion")
xlim([1010, 1090])
