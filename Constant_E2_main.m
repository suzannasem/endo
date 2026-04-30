%%%% Incorporating E2 as a constant input %%%%
%%% Simulation file %%%

%% Initial conditions %%

m00 = 9.0 * 10^5 ;
m10 = 100 ;
m20 = 4.5 * 10^4 ;
k00 = 5 * 10^5 ; 
ka0 = 2 * 10^4 ; 
e00 = 0 ; ef0 = 0; ea0 = 0;
Ra0 = 1; Rb0 = 0; 

init = [m00, m10, m20, k00, ka0, e00, ef0, ea0, Ra0, Rb0] ;

opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

%% Fitting E2 and P4 data

% grab data from table
McLKeefe = readtable('McL_Keefe.dat');

% assign columns to hormones
E2_data = McLKeefe.E2;
P4_data = McLKeefe.P4;
t = McLKeefe.time;

% "average" hormonal level
E2_mean = mean(E2_data) ;
P4_mean = mean(P4_data) ;

%% ODE Simulations %%

tspan = 0:0.5:1200 ;
[t,Y] = ode15s(@(t,y) Constant_E2_fn(t,y, 1e-6, 1e-5, E2_mean, P4_mean, 3e3, 4e3),tspan,init, opts); % disease free

% different values of K1, K2
[t1,Y1] = ode15s(@(t,y) Constant_E2_fn(t,y, 1e-6, 1e-5, E2_mean, P4_mean, 4e3, 3e3),tspan,init, opts); % slight malfunction
[t2,Y2] = ode15s(@(t,y) Constant_E2_fn(t,y, 1e-6, 1e-5, E2_mean, P4_mean, 5e3, 3e3),tspan,init, opts); % early disease
[t3,Y3] = ode15s(@(t,y) Constant_E2_fn(t,y, 1e-6, 1e-5, E2_mean, P4_mean, 6e3, 3e3),tspan,init, opts); % disease

% default (scenario 1)
m0 = Y(:,1);
m1 = Y(:,2);
m2 = Y(:,3);
k0 = Y(:,4);
ka = Y(:,5);
e0 = Y(:,6);
ef = Y(:,7);
ea = Y(:,8);
Ra = Y(:,9);
Rb = Y(:,10);

% early malfunction
m01 = Y1(:,1);
m11 = Y1(:,2);
m21 = Y1(:,3);
k01 = Y1(:,4);
ka1 = Y1(:,5);
e01 = Y1(:,6);
ef1 = Y1(:,7);
ea1 = Y1(:,8);
Ra1 = Y1(:,9);
Rb1 = Y1(:,10);

%  borderline
m02 = Y2(:,1);
m12 = Y2(:,2);
m22 = Y2(:,3);
k02 = Y2(:,4);
ka2 = Y2(:,5);
e02 = Y2(:,6);
ef2 = Y2(:,7);
ea2 = Y2(:,8);
Ra2 = Y2(:,9);
Rb2 = Y2(:,10);

% disease
m03 = Y3(:,1);
m13 = Y3(:,2);
m23 = Y3(:,3);
k03 = Y3(:,4);
ka3 = Y3(:,5);
e03 = Y3(:,6);
ef3 = Y3(:,7);
ea3 = Y3(:,8);
Ra3 = Y3(:,9);
Rb3 = Y3(:,10);


%% Plotting %%

% healthy profile

figure(101)
tc1 = tiledlayout(2,2);

nexttile(tc1)
plot(t, k0, '--', t, ka, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,400e3])
title("Natural Killer Cells")

nexttile(tc1)
plot(t, e0, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0, 15e3])
title("Eutopic Endometrial Cells")

nexttile(tc1)
semilogy(t, m0, '--', t, m1, '-', t, m2, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
title("Macrophage")

nexttile(tc1)
plot(t, ea, '-', t, ef, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells")

% diseased profilek1 = 6e3, k2 = 3e3
figure(106)
tc1 = tiledlayout(2,2) ;

nexttile(tc1)
plot(t3, k03, '--', t3, ka3, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,400e3])
title("Natural Killer Cells")

nexttile(tc1)
plot(t3, e03, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0, 15e3])
title("Eutopic Endometrial Cells")

nexttile(tc1)
semilogy(t3, m03, '--', t3, m13, '-', t3, m23, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
title("Macrophage")

nexttile(tc1)
plot(t3, ea3, '-', t3, ef3, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells")
%%
figure(203030303)
tc1 = tiledlayout(2,3) ;

nexttile(tc1)
semilogy(t, k0, '--', t, ka, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,400e3])
title("Natural Killer Cells in DF")

nexttile(tc1)
semilogy(t, m0, '--', t, m1, '-', t, m2, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel({})
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
title("Macrophage in DF")

nexttile(tc1)
plot(t, ea, '-', t, ef, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells in DF")

nexttile(tc1)
semilogy(t3, k03, '--', t3, ka3, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,400e3])
title("Natural Killer Cells in Disease")

nexttile(tc1)
semilogy(t3, m03, '--', t3, m13, '-', t3, m23, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel({})
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
title("Macrophage in Disease")

nexttile(tc1)
plot(t3, ea3, '-', t3, ef3, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells in Disease")

sgtitle('Immune and Uterine Cell Profiles in Disease & Disease-Free (DF) States', 'FontSize', 16)

%%
figure(202)
tc1 = tiledlayout(2,1);

nexttile(tc1)
plot(t, Ra, '--', 'LineWidth', 2, 'Color', '#5F2399', 'DisplayName', 'ER_\alpha')
hold on
plot(t, Rb, '-', 'LineWidth', 2, 'Color', '#3A5ABD', 'DisplayName', 'ER_\beta')
hold off
%xlabel("Time in days")
ylabel("Receptors / mL")
xlim([1006, 1056])
%xticks(1006:10:1073)
%xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Estrogen Receptors in DF")

nexttile(tc1)
line1 = plot(t3, Ra3, '--', 'LineWidth', 2, 'Color', '#5F2399', 'DisplayName', 'ER_\alpha');
hold on
line11 = plot(t3, Rb3, '-', 'LineWidth', 2, 'Color', '#3A5ABD', 'DisplayName', 'ER_\beta');
hold off
xlabel("Time in days")
ylabel("Receptors / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Estrogen Receptors in Disease")
fontsize(gcf, 12, "points")
hL = legend([line1, line11]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
hL.Title.String = 'Estrogen Receptor';
%title(tc1, 'Modeled and Actual Hormonal Levels During the Cycle', 'FontSize', 16)
%subtitle(tc1, "Disease-Free State")


%subplot(2,1,2)
%plot(t, Da, '-', t, Dab, '-', 'LineWidth', 2)
%legend("D_a", "D_{ab}")
%xlabel("Time in days")
%ylabel("Dimers / mL")
%xlim([1010, 1072])
%ylim([0,0.7e3])
%title("Estrogen Receptor Dimers on Macrophages")
%%
% ERb/ERa ratio (important in distinguishing disease states)

ratio = Rb ./ Ra ;
ratio1 = Rb1 ./ Ra1 ;
ratio2 = Rb2 ./ Ra2 ;
ratio3 = Rb3 ./ Ra3 ;

figure(203)
plot(t, Rb, '-', 'LineWidth',2)
hold on
plot(t1, Rb1, '--', 'LineWidth',2)
plot(t2, Rb2, ':', 'LineWidth',2)
plot(t3, Rb3, '-.', 'LineWidth',2)
hold off
legend("Disease-Free", "Slight Malfunction", "Early Disease", "Late Disease")
xlabel('Time in Days')
title('ER_\beta Amounts (Dimensionless)')
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,8e3])


% unbound ERb / ERa ratio
%erb_u = K1 - Dab - Rb ;
%era_u = K2 - Dab - Da - Ra ;

%figure(2000202092)
%plot(t, erb_u)
%hold on
%plot(t, era_u)
%% Different Disease States

figure(637)
tc1 = tiledlayout(2,1) ;

nexttile(tc1)
plot(t, Rb, '-', 'LineWidth', 2)
%hold on
%plot(t1, Rb1, '-', 'LineWidth', 2)
%hold off
%legend("1", "2")
xlabel("Time in days")
ylabel("Receptors / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.05])
title("ER_\beta Proportion for Disease Free States")

nexttile(tc1)
plot(t2, Rb2, '-', 'LineWidth', 2)
%hold on
%plot(t3, Rb3, '-', 'LineWidth', 2)
%hold off
%legend("3", "4")
xlabel("Time in days")
ylabel("Receptors / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0, 2001])
title("ER_\beta Proportion for Disease States")

figure(638)

subplot(2,1,1)
plot(t, ea, '-', 'LineWidth', 2)
hold on
plot(t1, ea1, '-', 'LineWidth', 2)
hold off
legend("1", "2")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.04])
title("Lesion Cell Count")

subplot(2,1,2)
plot(t2, ea2, '-', 'LineWidth', 2)
hold on
plot(t3, ea3, '-', 'LineWidth', 2)
hold off
legend("3", "4")
xlabel("Time in days")
ylabel("Receptors / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0, 1e9])

%% Genetic Tendency for ERb


% sequences of parameter values to test
K1_list = linspace(1e3,1e5,20);
K2_list = linspace(1e3,1e5,20);

% empty matrix for storing ea
ea_vals = zeros(20,20) ;

% empty matrix for storing erb
erb_vals = zeros(20,20) ;

% timing the simulation
tic

% loop to fill matrix
for j = 1:length(K1_list)
    for i = 1:length(K2_list)
            % run simulations
            [t_ea,Y_ea] = ode15s(@(t,y) Constant_E2_fn(t,y, 1e-6, 1e-5, E2_mean, P4_mean, K1_list(j), K2_list(i)),tspan,init, opts);
            % storing ea values for heatmap 1
            ea_vals(j,i) = Y_ea(end,8);
            % storing m0 values for proportions
            erb_vals(j,i) = Y_ea(end,10);
    end
end

% timing
elapsedTime = toc;  % stop timer
fprintf('Loop runtime: %.2f seconds\n', elapsedTime);
%%
% plot heatmap 
figure(4)

subplot(1,2,1)
imagesc([1e3,1e5],[1e3,1e5],log10(ea_vals))
set(gca, 'XScale', 'log', 'YScale', 'log','YDir','normal')
ylabel('K_1')
xlabel('K_2')
cb = colorbar ;
cb.Label.String = 'Lesion Cell Count in Powers of 10' ;
title("Lesion Cell Count")

subplot(1,2,2)
imagesc([1e3,1e5],[1e3,1e5],log10(erb_vals))
set(gca, 'XScale', 'log', 'YScale', 'log','YDir','normal')
ylabel('K_1')
xlabel('K_2')
cb = colorbar ;
cb.Label.String = 'ER_\beta Count, Dimensionless' ;
title('ER_\beta Count')

%% Different Amounts of Estrogen

% fix K2 for ratio calculation
K2_fixed = 1e4; 
ratio_list = logspace(-1,3,20); % 0.1 to 1000
E2_list = logspace(0,5,20);

% empty matrix for storing ea
ea_vals_E2 = zeros(length(E2_list), length(ratio_list));

% empty matrix for storing erb
erb_vals_E2 = zeros(length(E2_list), length(ratio_list));

% timing the simulation
tic

for r = 1:length(E2_list)
    
    for j = 1:length(ratio_list)
        
        ratio = ratio_list(j);
        K1 = ratio * K2_fixed;
        K2 = K2_fixed;
        
        [t_ea_E2,Y_ea_E2] = ode15s(@(t,y) Constant_E2_fn(t,y, ...
            1e-6, 1e-5, E2_list(r), P4_mean, K1, K2), ...
            tspan, init, opts);
        
        % Use final steady-state value
        ea_vals_E2(r,j) = Y_ea_E2(end,8);
        erb_vals_E2(r,j) = Y_ea_E2(end,10);
        
    end
end

% timing
elapsedTime = toc;  % stop timer
fprintf('Loop runtime: %.2f seconds\n', elapsedTime);

%%
% plot heatmap

figure(5)

subplot(1,2,1)
imagesc(ratio_list, E2_list,log10(ea_vals_E2))
set(gca, 'XScale', 'log', 'YScale', 'log','YDir','normal')
xlabel('K1 / K2 Ratio')
ylabel('E_2')
title('Lesion Cell Count')
cb = colorbar ;
cb.Label.String = 'Lesion Cell Count in Powers of 10' ;

subplot(1,2,2)
imagesc(ratio_list, E2_list, log10(erb_vals_E2))
set(gca, 'XScale', 'log', 'YScale', 'log','YDir','normal')
xlabel('K1 / K2 Ratio')
ylabel('E_2')
title('ER_\beta Count')
cb = colorbar ;
cb.Label.String = 'ER_\beta Count in Powers of 10' ;