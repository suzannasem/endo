%% Main FrankenFile %%

%% Fit data 

% grab data from table
McLKeefe = readtable('McL_Keefe.dat');

% assign columns to hormones
E2_data = McLKeefe.E2;
P4_data = McLKeefe.P4;
FSH_data = McLKeefe.FSH;
LH_data = McLKeefe.LH ;
time = McLKeefe.time + 1093;
% t = McLKeefe.time;

%% Initial Conditions

% initial conditions
FSHp0 = 117.86;
FSH0 = 154.6;
LHp0 = 251.73;
LH0 = 5.0888;
PH0 = 0.50465;
OM0 = 9.7718;
LB0 = 4.0704;
S0 = 1;
E20 = 41.047;
%E20 = 100 ;
P40 = 0.82942;
m00 = 9.0 * 10^5 ;
m10 = 100 ;
m20 = 4.5 * 10^4 ;
k00 = 5 * 10^5 ; 
ka0 = 2 * 10^4 ; 
e00 = 0 ; ef0 = 0; ea0 = 0;
Ra0 = 1;  Rb0 = 0; 

init1 = [m00, m10, m20, k00, ka0, e00, ef0, ea0, Ra0, Rb0, FSHp0, FSH0, LHp0, LH0, PH0, OM0, LB0, S0, E20, P40] ;
scale = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.0104, 0.97083, 0.12339, 1.0207, 0.074126, 0.87525, 0.54833, 0.11515, 0.23555, 1.636] ;

init = init1 .* scale ; % scaling

% timespan vector
tspan = 0:0.5:1210 ;

% options
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
%% Simulate

tic
% psi = 6e-8, 7e-8 gave a long cycle, 6e-7 gave numerical errors
% 9e-8 long cycle (~65 days?)
% 1e-7 - no full shedding (e0 doesn't hit 0) and long cycle length 
% 2e-7 - flat hormonal cycle, period gets smaller and disappears

% run
[t,Y] = ode23s(@(t,y) Franken_fn(t,y, 1e-6, 1e-5, 1e3, 4e3, 2e-7), tspan,init, opts); % healthy baseline
[t2,Y2] = ode23s(@(t,y) Franken_fn(t,y, 1e-6, 1e-5, 6e3, 3e3, 1e-10),tspan,init, opts); % disease + low E2 production by lesions
[t3,Y3] = ode23s(@(t,y) Franken_fn(t,y, 1e-6, 1e-5, 6e3, 3e3, 2e-7),tspan,init, opts); % disease + infertility
[t4,Y4] = ode23s(@(t,y) Franken_fn(t,y, 1e-6, 1e-5, 6e3, 3e3, 1e-7),tspan,init, opts); % disease + long cycle

% rename variables
% healthy
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
FSHp = Y(:,11);
FSH = Y(:,12);
LHp = Y(:,13);
LH = Y(:,14);
PH = Y(:,15);
OM = Y(:,16);
LB = Y(:,17);
S = Y(:,18);
E2 = Y(:,19);
P4 = Y(:,20);

% disease + low E2 prod by lesions
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
FSHp2 = Y2(:,11);
FSH2 = Y2(:,12);
LHp2 = Y2(:,13);
LH2 = Y2(:,14);
PH2 = Y2(:,15);
OM2 = Y2(:,16);
LB2 = Y2(:,17);
S2 = Y2(:,18);
E22 = Y2(:,19);
P42 = Y2(:,20);

% disease + high E2 prod
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
FSHp3 = Y3(:,11);
FSH3 = Y3(:,12);
LHp3 = Y3(:,13);
LH3 = Y3(:,14);
PH3 = Y3(:,15);
OM3 = Y3(:,16);
LB3 = Y3(:,17);
S3 = Y3(:,18);
E23 = Y3(:,19);
P43 = Y3(:,20);

% disease + medium E2 prod
m04 = Y4(:,1);
m14 = Y4(:,2);
m24 = Y4(:,3);
k04 = Y4(:,4);
ka4 = Y4(:,5);
e04 = Y4(:,6);
ef4 = Y4(:,7);
ea4 = Y4(:,8);
Ra4 = Y4(:,9);
Rb4 = Y4(:,10);
FSHp4 = Y4(:,11);
FSH4 = Y4(:,12);
LHp4 = Y4(:,13);
LH4 = Y4(:,14);
PH4 = Y4(:,15);
OM4 = Y4(:,16);
LB4 = Y4(:,17);
S4= Y4(:,18);
E24 = Y4(:,19);
P44 = Y4(:,20);

toc
%% Plotting

% Healthy system
figure(1)
tc1 = tiledlayout(2,2) ;

nexttile(tc1)
plot(t, k0, '--', t, ka, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,400e3])
title("Natural Killer Cells")

nexttile(tc1)
plot(t, e0, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0, 15e3])
title("Eutopic Endometrial Cells")

nexttile(tc1)
semilogy(t, m0, '--', t, m1, '-', t, m2, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
title("Macrophage")

nexttile(tc1)
plot(t, ea, '-', t, ef, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells")

title(tc1, "Immune and Endometrial Cell Profiles")
subtitle(tc1, "Disease-Free State")

figure(2);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t, FSH, '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Predicted');
line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Actual') ;
hold off
xticklabels({})
xlim([1100, 1153])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line2 = plot(t, LH, '-', 'Color', '#5F2399', 'LineWidth', 2);
line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1153])
title('LH (IU/L)');

nexttile(tcl)
hold on
line3 = plot(t, E2, '-', 'Color', '#5F2399', 'LineWidth', 2);
line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1153])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line4 = plot(t, P4, '-', 'Color', '#5F2399', 'LineWidth', 2);
line44 = plot(time, P4_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xlim([1100, 1153])
xticks(1100:10:1153)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
xlabel("Time (Days)")
title('P_4 (pg/mL)');

fontsize(gcf, 12, "points")
hL = legend([line1, line11]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
title(tcl, 'Modeled and Actual Hormonal Levels During the Cycle', 'FontSize', 16)
subtitle(tcl, "Disease-Free State")

figure(3)

subplot(2,2,1)
plot(t2, k02, '--', t2, ka2, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,400e3])
title("Natural Killer Cells")

subplot(2,2,2)
plot(t2, e02, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0, 15e3])
title("Eutopic Endometrial Cells")

subplot(2,2,3)
semilogy(t2, m02, '--', t2, m12, '-', t2, m22, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
title("Macrophage")

subplot(2,2,4)
plot(t2, ea2, '-', t2, ef2, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells")

figure(4);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t2, FSH2, '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Predicted');
line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Actual') ;
hold off
xticklabels({})
xlim([1100, 1183])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line2 = plot(t2, LH2, '-', 'Color', '#5F2399', 'LineWidth', 2);
line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1183])
title('LH (IU/L)');

nexttile(tcl)
hold on
line3 = plot(t2, E22, '-', 'Color', '#5F2399', 'LineWidth', 2);
line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1183])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line4 = plot(t2, P42, '-', 'Color', '#5F2399', 'LineWidth', 2);
line44 = plot(time, P4_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
xlabel("Time (Days)")
title('P_4 (pg/mL)');

fontsize(gcf, 12, "points")
hL = legend([line1, line11]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
sgtitle('Modeled and Actual Hormonal Levels During the Cycle, Mild Disease', 'FontSize', 16)

figure(10000)

subplot(2,2,1)
plot(t4, k04, '--', t4, ka4, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,400e3])
title("Natural Killer Cells")

subplot(2,2,2)
plot(t4, e04, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0, 15e3])
title("Eutopic Endometrial Cells")

subplot(2,2,3)
semilogy(t4, m04, '--', t4, m14, '-', t4, m24, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
title("Macrophage")

subplot(2,2,4)
plot(t4, ea4, '-', t4, ef4, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells")

time1 = time + 40 ;
figure(2000909);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t4, FSH4, '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Predicted');
%line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Actual') ;
hold off
xticklabels({})
xlim([1100, 1183])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line2 = plot(t4, LH4, '-', 'Color', '#5F2399', 'LineWidth', 2);
%line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1183])
title('LH (IU/L)');

nexttile(tcl)
hold on
line3 = plot(t4, E24, '-', 'Color', '#5F2399', 'LineWidth', 2);
%line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1183])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line4 = plot(t4, P44, '-', 'Color', '#5F2399', 'LineWidth', 2);
%line44 = plot(time, P4_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xlim([1100, 1183])
xticks(1100:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
xlabel("Time (Days)")
title('P_4 (pg/mL)');

fontsize(gcf, 12, "points")
hL = legend([line1, line11]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
sgtitle('Modeled and Actual Hormonal Levels During the (Long) Cycle', 'FontSize', 16)

figure(5)

subplot(2,2,1)
plot(t3, k03, '--', t3, ka3, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1140, 1200])
%ylim([0,400e3])
title("Natural Killer Cells")

subplot(2,2,2)
plot(t3, e03, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1140, 1200])
ylim([0, 1900])
title("Eutopic Endometrial Cells")

subplot(2,2,3)
semilogy(t3, m03, '--', t3, m13, '-', t3, m23, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1140, 1200])
title("Macrophage")

subplot(2,2,4)
plot(t3, ea3, '-', t3, ef3, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1140, 1200])
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells")

figure(6);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t3, FSH3, '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Predicted');
line11 = plot(time1, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Actual') ;
hold off
xticklabels({})
xlim([1140, 1200])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line2 = plot(t3, LH3, '-', 'Color', '#5F2399', 'LineWidth', 2);
line22 = plot(time1, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1140, 1200])
title('LH (IU/L)');

nexttile(tcl)
hold on
line3 = plot(t3, E23, '-', 'Color', '#5F2399', 'LineWidth', 2);
line33 = plot(time1, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1140, 1200])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line4 = plot(t3, P43, '-', 'Color', '#5F2399', 'LineWidth', 2);
line44 = plot(time1, P4_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xlim([1140, 1200])
xticks(1140:10:1200)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
xlabel("Time (Days)")
title('P_4 (pg/mL)');

fontsize(gcf, 12, "points")
hL = legend([line1, line11]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
sgtitle('Modeled and Actual Hormonal Levels in Disease State', 'FontSize', 16)

figure(7)
subplot(2,1,1)
plot(t, Ra, '-', t, Rb, '-', 'LineWidth', 2)
legend("ER_\alpha", "ER_\beta")
xlabel("Time in days")
ylabel("Receptors / mL")
xlim([1010, 1072])
%ylim([0,0.7e3])
title("Estrogen Receptors on Macrophages")

subplot(2,1,2)
plot(t3, Ra3, '-', t3, Rb3, '-', 'LineWidth', 2)
legend("ER_\alpha", "ER_\beta")
xlabel("Time in days")
ylabel("Receptors / mL")
xlim([1010, 1072])
%ylim([0,0.7e3])
title("Estrogen Receptors on Macrophages in Disease")


ratio = Rb ./ Ra ;
ratio3 = Rb3 ./ Ra3 ;

figure(8)
plot(t, ratio, '-', 'LineWidth',2)
hold on
plot(t3, ratio3, '-.', 'LineWidth',2)
hold off
legend("Disease-Free", "Late Disease")
xlabel('Time in Days')
title('ER_\beta to ER_\alpha Ratio')
xlim([1010, 1072])
%%
figure(1110)
plot(t3, e03, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([200, 1200])
xticks(200:100:1200)
xticklabels({"0", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1100", "1200"})
%ylim([0, 15e3])
title("Eutopic Endometrial Cells in Sterile State")
subtitle("High E_2 synthesis by Lesions")

%%
figure(9);
tcl = tiledlayout(3, 1);

nexttile(tcl)
hold on
%line1 = plot(t, ph_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t, PH, '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Disease-Free');
%line13 = plot(t, ph_vals(5,:), '-','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t3, PH3, '-','Color', '#480082','LineWidth',2, 'DisplayName','Diseased');
%xticklabels({})
xlim([1090, 1183])
title('Follicular Phase');

nexttile(tcl)
hold on
%line1 = plot(t, om_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t, OM, '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Disease-Free');
%line13 = plot(t, om_vals(5,:), '-','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t3, OM3, '-','Color', '#480082','LineWidth',2, 'DisplayName','Disease');
%xticklabels({})
xlim([1090, 1183])
title('Ovulatory Phase');

nexttile(tcl)
hold on
%line1 = plot(t, lb_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t, LB, '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Disease-Free');
%line13 = plot(t, lb_vals(5,:), '-','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t3, LB3, '-','Color', '#480082','LineWidth',2, 'DisplayName','Disease');
%xticklabels({})
xlim([1090, 1183])
title('Luteal Phase');
xticks(1100:10:1160)
%xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})

fontsize(gcf, 12, "points")
hL = legend([line12, line14]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
sgtitle('Ovulatory Cycle', 'FontSize', 16)

%%

figure(2030303039)
tc1 = tiledlayout(3,3,'TileSpacing','compact') ;

nexttile(tc1)
semilogy(t2, k02, '--', t2, ka2, '-', 'LineWidth', 2)
legend("K_0", "K_A")
% xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1006, 1056])
% xticks(1006:10:1073)
xticklabels({})
%ylim([0,400e3])
% title("Natural Killer Cells in Mild Disease")

nexttile(tc1)
semilogy(t2, m02, '--', t2, m12, '-', t2, m22, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
% xlabel("Time in days")
ylabel({})
xlim([1006, 1056])
% xticks(1006:10:1073)
xticklabels({})
% title("Macrophage in Mild Disease")

nexttile(tc1)
semilogy(t2, ea2, '-', t2, ef2, '-', 'LineWidth', 2)
legend("E_A", "E_F")
% xlabel("Time in days")
% ylabel("Cells / mL")
xlim([1006, 1056])
% xticks(1006:10:1073)
xticklabels({})
%ylim([0,0.7e3])
% title("Ectopic Endometrial Cells Mild Disease")

nexttile(tc1)
semilogy(t4, k04, '--', t4, ka4, '-', 'LineWidth', 2)
legend("K_0", "K_A")
% xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1006, 1056])
% xticks(1006:10:1073)
xticklabels({})
%ylim([0,400e3])
% title("Natural Killer Cells in Moderate Disease")

nexttile(tc1)
semilogy(t4, m04, '--', t4, m14, '-', t4, m24, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
% xlabel("Time in days")
ylabel({})
xlim([1006, 1056])
% xticks(1006:10:1073)
xticklabels({})
% title("Macrophage in Moderate Disease")

nexttile(tc1)
semilogy(t4, ea4, '-', t4, ef4, '-', 'LineWidth', 2)
legend("E_A", "E_F")
% xlabel("Time in days")
% ylabel("Cells / mL")
xlim([1006, 1056])
% xticks(1006:10:1073)
xticklabels({})
%ylim([0,0.7e3])
% title("Ectopic Endometrial Cells in Moderate Disease")

nexttile(tc1)
semilogy(t3, k03, '--', t3, ka3, '-', 'LineWidth', 2)
legend("K_0", "K_A")
xlabel("Time in days", 'FontSize', 16)
ylabel("Cells / mL (log scale)")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,400e3])
% title("Natural Killer Cells in Sterile Disease")

nexttile(tc1)
semilogy(t3, m03, '--', t3, m13, '-', t3, m23, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days", 'FontSize', 16)
ylabel({})
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
% title("Macrophage in Sterile Disease")

nexttile(tc1)
semilogy(t3, ea3, '-', t3, ef3, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days", 'FontSize', 16)
% ylabel("Cells / mL")
xlim([1006, 1056])
xticks(1006:10:1073)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
%ylim([0,0.7e3])
% title("Ectopic Endometrial Cells in Sterile Disease")

% Row Titles 
annotation('textbox', [.01 .75 .05 .1], 'String', 'Fertile Disease', 'EdgeColor', 'none', 'FontAngle', 'italic', 'FontSize', 15);
annotation('textbox', [.01 .47 .05 .1], 'String', 'Long Cycle Disease', 'EdgeColor', 'none', 'FontAngle', 'italic', 'FontSize', 15);
annotation('textbox', [.01 .15 .05 .1], 'String', 'Sterile Disease', 'EdgeColor', 'none', 'FontAngle', 'italic', 'FontSize', 15);

% Column Titles
annotation('textbox', [.2 .9 .1 .05], 'String', 'NK Cells', 'EdgeColor', 'none', 'FontAngle', 'italic', 'FontSize', 15);
annotation('textbox', [.47 .9 .1 .05], 'String', 'Macrophage', 'EdgeColor', 'none', 'FontAngle', 'italic', 'FontSize', 15);
annotation('textbox', [.71 .9 .2 .05], 'String', 'Ectopic Endometrial Cells', 'EdgeColor', 'none', 'FontAngle', 'italic', 'FontSize', 15);

% titles and formatting
% set(gca, 'FontWeight', 'bold', 'FontSize', 14)
sgtitle('Immune and Endometrial Cell Profiles in Disease States', 'FontSize', 22);

%% Amount of E2 produced by Lesions
% sequences of parameter values to test 
psi_list = linspace(1e-10,1.5e-7,5);

tic

for j = 1:length(psi_list)

    [t,Y] = ode23s(@(t,y) Franken_fn(t,y, 1e-6, 1e-5, 6e3, 3e3, psi_list(j)), ...
                   tspan, init, opts);

    % [t,Y] = ode23s(@(t,y) Franken_fn(t,y, 1e-6, 1e-5, 6e3, 3e3, psi_list(j)), tspan,init, opts);

    % Preallocate after first run
    if j == 1
        Nt = length(t);

        fshp_vals = zeros(length(psi_list), Nt);
        fsh_vals  = zeros(length(psi_list), Nt);
        lhp_vals  = zeros(length(psi_list), Nt);
        lh_vals   = zeros(length(psi_list), Nt);
        ph_vals   = zeros(length(psi_list), Nt);
        om_vals   = zeros(length(psi_list), Nt);
        lb_vals   = zeros(length(psi_list), Nt);
        e2_vals   = zeros(length(psi_list), Nt);
        p4_vals   = zeros(length(psi_list), Nt);
        e0_vals = zeros(length(psi_list), Nt);
       % fprintf('psi = %.2e, max FSH = %.2e\n', psi_list(j), max(Y(:,2)));
    end

    % disp(init)
    % Store full trajectories (transpose to make row)
    fshp_vals(j,:) = Y(:,11)';
    fsh_vals(j,:)  = Y(:,12)';
    lhp_vals(j,:)  = Y(:,13)';
    lh_vals(j,:)   = Y(:,14)';
    ph_vals(j,:)   = Y(:,15)';
    om_vals(j,:)   = Y(:,16)';
    lb_vals(j,:)   = Y(:,17)';
    e2_vals(j,:)   = Y(:,19)';
    p4_vals(j,:)   = Y(:,20)';
    e0_vals(j,:)   = Y(:,6)';
end

toc

% timing
elapsedTime = toc;  % stop timer
% fprintf('Loop runtime: %.2f seconds\n', elapsedTime);

%%

figure(10);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t, fsh_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t, fsh_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, fsh_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line14 = plot(t, fsh_vals(10,:), '-','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Data') ;
hold off
xticklabels({})
xlim([1100, 1200])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line2 = plot(t, lh_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line22 = plot(t, lh_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line23 = plot(t, lh_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line24 = plot(t, lh_vals(10,:), '-','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xticklabels({})
xlim([1100, 1200])
title('LH (IU/L)');

nexttile(tcl)
hold on
line1 = plot(t, e2_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t, e2_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, e2_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line14 = plot(t, e2_vals(10,:), '-','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xticklabels({})
xlim([1100, 1200])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line1 = plot(t, p4_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t, p4_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, p4_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line14 = plot(t, p4_vals(10,:), '-','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line44 = plot(time, P4_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xlim([1100, 1200])
xticks(1100:10:1200)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90"})
xlabel("Time (Days)")
title('P_4 (pg/mL)');

%hL = legend([line1, line12, line13]);
%hL.FontSize = 14;
%hL.Layout.Tile = 'east';
%hL.Title.String = 'E_2 Production Rate';
fontsize(gcf, 12, "points")

title(tcl, 'Hormonal Levels throughout the Cycle', 'FontSize', 16);
subtitle(tcl, 'Varying E_2 Production Rate by Lesions');
%%
figure(11)
tcl = tiledlayout(3, 1);

nexttile(tcl)
hold on
line1 = plot(t, ph_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t, ph_vals(3,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, ph_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line14 = plot(t, fsh_vals(10,:), '-','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Data') ;
hold off
xticklabels({})
xlim([1100, 1200])
title('Follicular Phase');

nexttile(tcl)
hold on
line2 = plot(t, om_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line22 = plot(t, om_vals(3,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line23 = plot(t, om_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line24 = plot(t, lh_vals(10,:), '-','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xticklabels({})
xlim([1100, 1200])
title('Ovulatory Phase');

nexttile(tcl)
hold on
line1 = plot(t, lb_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t, lb_vals(3,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, lb_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line14 = plot(t, e2_vals(10,:), '-','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xticklabels({})
xlim([1100, 1200])
title('Luteal Phase');
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90"})
xlabel("Time (Days)")

hL = legend([line1, line12, line13]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
hL.Title.String = 'E_2 Production Rate';
fontsize(gcf, 12, "points")

title(tcl, 'Ovulatory Cycle', 'FontSize', 16);
subtitle(tcl, 'Varying E_2 Production Rate by Lesions');

%% 
figure(22222)
plot(t, e0_vals(3,:), 'LineWidth', 2, 'Color', '#BF40BF')
xticklabels({})
xlim([1000, 1200])
title('Eutopic Endometrial Cells with Moderate E_2 Synthesis');
xticks(1000:30:1200)
xticklabels({"0", "30", "60", "90", "120", "150", "180", "210", "240", "270", "300"})
xlabel("Time (Days)")