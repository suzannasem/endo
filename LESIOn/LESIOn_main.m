%% Main LESIOn file %%

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
LE20 = 0;

init1 = [m00, m10, m20, k00, ka0, e00, ef0, ea0, Ra0, Rb0, FSHp0, FSH0, LHp0, LH0, PH0, OM0, LB0, S0, E20, P40, LE20] ;
scale = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.0104, 0.97083, 0.12339, 1.0207, 0.074126, 0.87525, 0.54833, 0.11515, 0.23555, 1.636, 1] ;

init = init1 .* scale ; % scaling

% timespan vector
tspan = 0:0.5:1210 ;

% parameters
p = LESIOn_params();

% options
opts = odeset('RelTol',1e-6,'AbsTol',1e-9); 
%% Simulate

tic

% run
[t,Y] = ode15s(@(t,y) LESIOn_fn(t,y, 1e3, 4e3, 2e-7, p), tspan,init, opts); % healthy baseline
[t2,Y2] = ode15s(@(t,y) LESIOn_fn(t,y, 6e3, 3e3, 1e-10, p),tspan,init, opts); % disease + low E2 production by lesions
[t3,Y3] = ode15s(@(t,y) LESIOn_fn(t,y, 6e3, 3e3, 1e-7, p),tspan,init, opts); % disease + high psi
[t4,Y4] = ode15s(@(t,y) LESIOn_fn(t,y, 6e3, 3e3, 1e-8, p),tspan,init, opts); % disease + moderate psi

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
LE2 = Y(:,21);

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
LE22 = Y2(:,21);

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
LE3 = Y3(:,21);

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
LE4 = Y4(:,21);

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
%%
figure(200);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t, FSH, '-','Color', '#F99F9F','LineWidth',2, 'DisplayName','Disease-Free');
line2 = plot(t2, FSH2, '--','Color', '#6A627F','LineWidth',2, 'DisplayName','Disease + Low');
%line5 = plot(t5, FSH5, '-.','Color', '#428576','LineWidth',2, 'DisplayName','Disease + Mid');
line4 = plot(t4, FSH4, ':','Color', '#86892A','LineWidth',2, 'DisplayName','Disease + High');
line3 = plot(t3, FSH3, '.','Color', '#5F2399','LineWidth',2, 'DisplayName','Disease + Very High');
%line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Actual') ;
hold off
xticklabels({})
xlim([1100, 1153])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line1 = plot(t, LH, '-','Color', '#F99F9F','LineWidth',2, 'DisplayName','Disease-Free');
line2 = plot(t2, LH2, '--','Color', '#6A627F','LineWidth',2, 'DisplayName','Disease + Low');
%line5 = plot(t5, LH5, '-.','Color', '#428576','LineWidth',2, 'DisplayName','Disease + Mid');
line4 = plot(t4, LH4, ':','Color', '#86892A','LineWidth',2, 'DisplayName','Disease + High');
line3 = plot(t3, LH3, '.','Color', '#5F2399','LineWidth',2, 'DisplayName','Disease + Very High');
%line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1153])
title('LH (IU/L)');

nexttile(tcl)
hold on
line1 = plot(t, E2, '-','Color', '#F99F9F','LineWidth',2, 'DisplayName','Disease-Free');
line2 = plot(t2, E22, '--','Color', '#6A627F','LineWidth',2, 'DisplayName','Disease + Low');
%line5 = plot(t5, E25, '-.','Color', '#428576','LineWidth',2, 'DisplayName','Disease + Mid');
line4 = plot(t4, E24, ':','Color', '#86892A','LineWidth',2, 'DisplayName','Disease + High');
line3 = plot(t3, E23, '.','Color', '#5F2399','LineWidth',2, 'DisplayName','Disease + Very High');
%line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1153])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line1 = plot(t, P4, '-','Color', '#F99F9F','LineWidth',2, 'DisplayName','Disease-Free');
line2 = plot(t2, P42, '--','Color', '#6A627F','LineWidth',2, 'DisplayName','Disease + Low');
%line5 = plot(t5, P45, '-.','Color', '#428576','LineWidth',2, 'DisplayName','Disease + Mid');
line4 = plot(t4, P44, ':','Color', '#86892A','LineWidth',2, 'DisplayName','Disease + High');
line3 = plot(t3, P43, '.','Color', '#5F2399','LineWidth',2, 'DisplayName','Disease + Very High');
%line44 = plot(time, P4_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xlim([1100, 1153])
xticks(1100:10:1153)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
xlabel("Time (Days)")
title('P_4 (pg/mL)');

fontsize(gcf, 12, "points")
hL = legend([line1, line2, line4, line3]);
hL.FontSize = 14;
%hL.Title = "Estradiol Synthesis Rate" ;
hL.Layout.Tile = 'east';
title(tcl, 'Hormonal Levels During the Cycle', 'FontSize', 16)
subtitle(tcl, "Varying \psi")
%%

function rgb = hex2rgb(hex)
    hex = erase(hex, '#');
    rgb = reshape(sscanf(hex,'%2x'),1,3) / 255;
end

colors = [
    hex2rgb('#E11584');
    hex2rgb('#8D4585');
    hex2rgb('#814C41');
    hex2rgb('#6F2DA8');
    hex2rgb('#9B43C7');
];

figure(201); clf
tcl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

lw = 2;
xlimits = [1100, 1153];

% ---------- Follicular ----------
nexttile
hold on
l(1) = plot(t,  PH,  '-', 'Color', colors(1,:), 'LineWidth', lw);
l(2) = plot(t2+1, PH2, '--', 'Color', colors(2,:), 'LineWidth', lw);
l(3) = plot(t4, PH4, '-.', 'Color', colors(4,:), 'LineWidth', lw);
l(4) = plot(t3, PH3, ':', 'Color', colors(5,:), 'LineWidth', lw);
hold off
title('Follicular Phase')
xlim(xlimits)
xticklabels([])
box off

% ---------- Ovulatory ----------
nexttile
hold on
plot(t,  OM,  '-', 'Color', colors(1,:), 'LineWidth', lw);
plot(t2+1, OM2, '--', 'Color', colors(2,:), 'LineWidth', lw);
plot(t4, OM4, '-.', 'Color', colors(4,:), 'LineWidth', lw);
plot(t3, OM3, ':', 'Color', colors(5,:), 'LineWidth', lw);
hold off
title('Ovulatory Phase')
xlim(xlimits)
xticklabels([])
box off

% ---------- Luteal ----------
nexttile
hold on
plot(t,  LB,  '-', 'Color', colors(1,:), 'LineWidth', lw);
plot(t2+1, LB2, '--', 'Color', colors(2,:), 'LineWidth', lw);
plot(t4, LB4, '-.', 'Color', colors(4,:), 'LineWidth', lw);
plot(t3, LB3, ':', 'Color', colors(5,:), 'LineWidth', lw);
hold off
title('Luteal Phase')
xlim(xlimits)
xticks(1100:10:1153)
xticklabels({"0","10","20","30","40","50","60","70"})
xlabel("Time (Days)")
box off
    
title(tcl, 'Ovarian Follicle Growth During the Cycle', 'FontSize', 16)
subtitle(tcl, "Varying E_2 Synthesis Rate by Lesions (\psi)")

% - - - - - Hormonal Cycle - - - - - - %
figure(200); clf
tcl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

lw = 2;
xlimits = [1100, 1153];

% ---------- FSH ----------
nexttile
hold on
l(1) = plot(t,  FSH,  '-', 'Color', colors(1,:), 'LineWidth', lw);
l(2) = plot(t2+1, FSH2, '--', 'Color', colors(2,:), 'LineWidth', lw);
l(3) = plot(t4, FSH4, '-.', 'Color', colors(4,:), 'LineWidth', lw);
l(4) = plot(t3, FSH3, ':', 'Color', colors(5,:), 'LineWidth', lw);
hold off
title('FSH (IU/L)')
xlim(xlimits)
xticklabels([])
box off

% ---------- LH ----------
nexttile
hold on
plot(t,  LH,  '-', 'Color', colors(1,:), 'LineWidth', lw);
plot(t2+1, LH2, '--', 'Color', colors(2,:), 'LineWidth', lw);
plot(t4, LH4, '-.', 'Color', colors(4,:), 'LineWidth', lw);
plot(t3, LH3, ':', 'Color', colors(5,:), 'LineWidth', lw);
hold off
title('LH (IU/L)')
xlim(xlimits)
xticklabels([])
box off

% ---------- E2 ----------
nexttile
hold on
plot(t,  E2,  '-', 'Color', colors(1,:), 'LineWidth', lw);
plot(t2+1, E22, '--', 'Color', colors(2,:), 'LineWidth', lw);
plot(t4, E24, '-.', 'Color', colors(4,:), 'LineWidth', lw);
plot(t3, E23, ':', 'Color', colors(5,:), 'LineWidth', lw);
hold off
title('E_2 (pg/mL)')
xlim(xlimits)
xticklabels([])
box off

% ---------- P4 ----------
nexttile
hold on
plot(t,  P4,  '-', 'Color', colors(1,:), 'LineWidth', lw);
plot(t2+1, P42, '--', 'Color', colors(2,:), 'LineWidth', lw);
plot(t4, P44, '-.', 'Color', colors(4,:), 'LineWidth', lw);
plot(t3, P43, ':', 'Color', colors(5,:), 'LineWidth', lw);
hold off
title('P_4 (pg/mL)')
xlim(xlimits)
xticks(1100:10:1153)
xticklabels({"0","10","20","30","40","50","60","70"})
xlabel("Time (Days)")
box off

% ---------- Legend ----------
hL = legend(l, ...
    {'Disease-Free','Disease: \psi = 1e-10','Disease: \psi = 1e-8','Disease: \psi = 1e-7'}, ...
    'FontSize',12);
hL.Layout.Tile = 'east';
hL.Title.String = 'Disease State';

title(tcl, 'Hormonal Levels During the Cycle', 'FontSize', 16)
subtitle(tcl, "Varying E_2 Synthesis Rate by Lesions (\psi)")

%%
% mild disease (no hormonal effects)
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
%%
% barely-there disease cycle
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

figure(41);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t4, FSH4, '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Predicted');
line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Actual') ;
hold off
xticklabels({})
xlim([1100, 1183])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line2 = plot(t4, LH4, '-', 'Color', '#5F2399', 'LineWidth', 2);
line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1183])
title('LH (IU/L)');

nexttile(tcl)
hold on
line3 = plot(t4, E24, '-', 'Color', '#5F2399', 'LineWidth', 2);
line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2);
hold off
xticklabels({})
xlim([1100, 1183])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line4 = plot(t4, P44, '-', 'Color', '#5F2399', 'LineWidth', 2);
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
%% 

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
%%
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
xlim([1000, 1200])
%ylim([0, 1900])
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
%%
% 
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


%% ER Calculations and Plots

calphaP = p.calphaP ;
calphaM = p.calphaM ;
cbetaP = p.cbetaP ;
cbetaM = p.cbetaM ;

% df: K1 = 1e3, K2 = 4e3
Da = p.calphaP / p.calphaM .* Ra.^2 * 4e3;
Dab = p.cbetaP / p.cbetaM .* 4e3 .* Ra .* Rb ;

% disease: K1 = 6e3, K2 = 3e3
Da2 = p.calphaP / p.calphaM .* Ra2.^2 * 3e3;
Dab2 = p.cbetaP / p.cbetaM .* 3e3 .* Ra2 .* Rb2 ;

Da3 = p.calphaP / p.calphaM .* Ra3.^2 * 3e3;
Dab3 = p.cbetaP / p.cbetaM .* 3e3 .* Ra3 .* Rb3 ;

Da4 = p.calphaP / p.calphaM .* Ra4.^2 * 3e3;
Dab4 = p.cbetaP / p.cbetaM .* 3e3 .* Ra4 .* Rb4 ;

% unbound receptors: Rb^u = K1 - Dab - Rb^b, Ra^u = K2 - Dab - Da - Ra^b

Rbu = 1 - Dab - Rb ;
Rau = 1 - (1e3/4e3) .* Dab - Da - Ra ;

Rbu2 = 1 - Dab2 - Rb2 ;
Rau2 = 1 - (6e3/3e3) .* Dab2 - Da2 - Ra2;

Rbu3 = 1 - Dab3 - Rb3 ;
Rau3 = 1 - (6e3/3e3) .*Dab3 - Da3 - Ra3;

Rbu4 = 1 - Dab4 - Rb4 ;
Rau4 = 1 - (6e3/3e3) .*Dab4 - Da4 - Ra4;

figure(35)
tc3 = tiledlayout(2,2);

nexttile(tc3)
hold on
line1 = plot(t, Rau+Ra, '-','Color', '#5F2399','LineWidth',2, 'DisplayName', 'ER_\alpha Ratio');
line11 = plot(t, Rbu+Rb, '--', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'ER_\beta Ratio') ;
hold off
xlim([1100, 1200])
ylim([0, 1e-5])
xticks({})
title("Disease-Free")

nexttile(tc3)
hold on
line1 = plot(t, Rau2+Ra2, '-','Color', '#5F2399','LineWidth',2, 'DisplayName', 'ER_\alpha Ratio');
line11 = plot(t, Rbu2+Rb2, '--', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'ER_\beta Ratio') ;
hold off
xlim([1100, 1200])
xticks({})
title("Disease: \psi = 1e-10")

nexttile(tc3)
hold on
line1 = plot(t, Rau4+Ra4, '-','Color', '#5F2399','LineWidth',2, 'DisplayName', 'ER_\alpha Ratio');
line11 = plot(t, Rbu4+Rb4, '--', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'ER_\beta Ratio') ;
hold off
xlim([1100, 1200])
xticks(1100:10:1200)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})
xlabel("Time (Days)")
title("Disease: \psi = 1e-8")

nexttile(tc3)
hold on
line1 = plot(t, Rau3+Ra3, '-','Color', '#5F2399','LineWidth',2, 'DisplayName', 'ER_\alpha');
line11 = plot(t, Rbu3+Rb3, '--', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'ER_\beta') ;
hold off
xlim([1100, 1200])
xticks(1100:10:1200)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})
xlabel("Time (Days)")
title("Disease: \psi = 1e-7")

fontsize(gcf, 12, "points")
hL = legend([line1, line11]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
sgtitle('ER_\alpha and ER_\beta Amounts (Dimensionless)', 'FontSize', 16)

% ratios?
rat1 = (Rau + Ra) ./ (Rbu + Rb) ;
rat2 = (Rau2 + Ra2) ./ (Rbu2 + Rb2) ;
rat3 = (Rau3 + Ra3) ./ (Rbu3 + Rb3) ;
rat4 = (Rau4 + Ra4) ./ (Rbu4 + Rb4) ;

%%
figure(44)

subplot(1,2,1)
plot(t, rat1, 'LineWidth', 2)
xlim([1100, 1200])
ylim([0, 0.5])
title("ER_\alpha to ER_\beta Ratio in DF")

subplot(1,2,2)
plot(t2, rat2,'-',  t4, rat4,':', t3, rat3,'--', 'LineWidth', 2)
legend("Low E_2 Production", "Moderate E_2 Production", "High E_2 Production")
xlim([1100, 1200])
ylim([0, 1e-5])
title("ER_\alpha to ER_\beta Ratio in Disease States")
%%
figure(1110)
plot(t3, e03, 'LineWidth', 2)
xlabel("Time in days")
ylabel("Cells / mL")
xlim([200, 1200])
xticks(200:100:1200)
xticklabels({"0", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1100", "1200"})
%ylim([0, 1])
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
line14 = plot(t3, PH3, '--','Color', '#480082','LineWidth',2, 'DisplayName','Sterile Disease');
line15 = plot(t4, PH4, '-.','Color', '#480082','LineWidth',2, 'DisplayName','Barely-Cycle Disease');
line16 = plot(t5, PH5, ':','Color', '#480082','LineWidth',2, 'DisplayName','Long-Cycle Disease');
%xticklabels({})
xlim([1090, 1183])
title('Follicular Phase');

nexttile(tcl)
hold on
%line1 = plot(t, om_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t, OM, '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Disease-Free');
%line13 = plot(t, om_vals(5,:), '-','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t3, OM3, '--','Color', '#480082','LineWidth',2, 'DisplayName','Sterile Disease');
line15 = plot(t4, OM4, '-.','Color', '#480082','LineWidth',2, 'DisplayName','Barely-Cycle Disease');
line16 = plot(t5, OM5, ':','Color', '#480082','LineWidth',2, 'DisplayName','Long-Cycle Disease');
%xticklabels({})
xlim([1090, 1183])
title('Ovulatory Phase');

nexttile(tcl)
hold on
%line1 = plot(t, lb_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t, LB, '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Disease-Free');
%line13 = plot(t, lb_vals(5,:), '-','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t3, LB3, '--','Color', '#480082','LineWidth',2, 'DisplayName','Sterile Disease');
line15 = plot(t4, LB4, '-.','Color', '#480082','LineWidth',2, 'DisplayName','Barely-Cycle Disease');
line16 = plot(t5, LB5, ':','Color', '#480082','LineWidth',2, 'DisplayName','Long-Cycle Disease');
%xticklabels({})
xlim([1090, 1183])
title('Luteal Phase');
xticks(1100:10:1160)
%xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})

fontsize(gcf, 12, "points")
hL = legend([line12, line14, line15, line16]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
sgtitle('Ovulatory Cycle', 'FontSize', 16)
%%
figure(1010101001)
tcg = tiledlayout(2,1);

% --- low amounts ---
nexttile(tcg)
plot(t, LE2, '-', 'Color', colors(1,:), 'LineWidth', 2)
hold on
plot(t2, LE2, '--', 'Color', colors(2,:), 'LineWidth', 2)
hold off

legend('Disease-Free', 'Disease: \psi = 1e-10', 'Location', 'best')
xlim([1090, 1183])
xticks([])
ylim([0, 13e-8])

% --- high amounts ---
nexttile(tcg)
plot(t4, LE4, '-.', 'Color', colors(4,:), 'LineWidth', 2)
hold on
plot(t3, LE3, ':', 'Color', colors(5,:), 'LineWidth', 2)
hold off

legend('Disease: \psi = 1e-8','Disease: \psi = 1e-7', 'Location', 'best')
xlim([1090, 1183])
xlabel('Time in Days')
xticks(1090:10:1183)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})

title(tcg, 'Local Estradiol Levels (pg/mL)')
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
psi_list = linspace(1e-10,1e-7,5);

tic

for j = 1:length(psi_list)

    [t,Y] = ode15s(@(t,y) LESIOn_fn(t,y, 6e3, 3e3, psi_list(j), p), ...
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

end

toc

% timing
elapsedTime = toc;  % stop timer
%fprintf('Loop runtime: %.2f seconds\n', elapsedTime);

%%

% note: changing psi creates a horizontal shift
% to make hormonal levels easier to compare, we shift some results
% horizontally (t+5)

figure(10);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t+5, fsh_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t+5, fsh_vals(2,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, fsh_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Data') ;
hold off
xticklabels({})
xlim([1100, 1200])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line2 = plot(t+5, lh_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line22 = plot(t+5, lh_vals(2,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line23 = plot(t, lh_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xticklabels({})
xlim([1100, 1200])
title('LH (IU/L)');

nexttile(tcl)
hold on
line1 = plot(t+5, e2_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t+5, e2_vals(2,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, e2_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
%line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xticklabels({})
xlim([1100, 1200])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line1 = plot(t+5, p4_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t+5, p4_vals(2,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, p4_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
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
line1 = plot(t+5, ph_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t+5, ph_vals(3,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t, ph_vals(4,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
line24 = plot(t-10, ph_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','Very High');
hold off
xticklabels({})
xlim([1100, 1200]) 
title('Follicular Phase');

nexttile(tcl)
hold on
line2 = plot(t+5, om_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line22 = plot(t+5, om_vals(3,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line23 = plot(t-10, om_vals(4,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
line24 = plot(t-10, om_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','Very High');
hold off
xticklabels({})
xlim([1100, 1200])
title('Ovulatory Phase');

nexttile(tcl)
hold on
line1 = plot(t+5, lb_vals(1,:), '--','Color', '#5F2399','LineWidth',2, 'DisplayName','Low');
line12 = plot(t+5, lb_vals(3,:), '-','Color', '#BF40BF','LineWidth',2, 'DisplayName','Moderate');
line13 = plot(t-10, lb_vals(4,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High');
line24 = plot(t-10, lb_vals(5,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','Very High');
hold off
xticklabels({})
xlim([1100, 1200])
title('Luteal Phase');
xticks(1100:10:1200)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90"})
xlabel("Time (Days)")

hL = legend([line1, line12, line13, line24]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
hL.Title.String = 'E_2 Production Rate';
fontsize(gcf, 12, "points")

title(tcl, 'Ovulatory Cycle', 'FontSize', 16);
subtitle(tcl, 'Varying E_2 Production Rate by Lesions');
