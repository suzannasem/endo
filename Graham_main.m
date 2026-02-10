%% Simulation File for Graham reduced model %%

%% Setup %%

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
P40 = 0.82942;

init = [FSHp0, FSH0, LHp0, LH0, PH0, OM0, LB0, S0, E20, P40] ;
scale = [1.0104, 0.97083, 0.12339, 1.0207, 0.074126, 0.87525, 0.54833, 0.11515, 0.23555, 1.636] ;

init = init .* scale ; % scaling

% timespan vector
tspan = 0:0.5:1210 ;

%% Simulation %%
[t,Y] = ode15s(@(t,y) Graham(t,y,graham_params), tspan,init); 

% rename variable
FSHp = Y(:,1);
FSH = Y(:,2);
LHp = Y(:,3);
LH = Y(:,4);
PH = Y(:,5);
OM = Y(:,6);
LB = Y(:,7);
S = Y(:,8);
E2 = Y(:,9);
P4 = Y(:,10);

%% Data from McKeefe

% grab data from table
McLKeefe = readtable('McL_Keefe.dat');

% assign columns to hormones
E2_data = McLKeefe.E2;
P4_data = McLKeefe.P4;
FSH_data = McLKeefe.FSH;
LH_data = McLKeefe.LH ;
time = McLKeefe.time + 1093; % shifting data points to match with model
% note: model is plotted from 1100 timesteps to capture periodic behavior
% after transient behavior


%% Plotting

figure(1)

subplot(4,1,1)
plot(t, FSH, '-',time, FSH_data,'x','LineWidth',2)
legend("Predicted", "Actual")
title("FSH")
xlim([1100, 1153])
xticklabels({})

subplot(4,1,2)
plot(t,LH, '--', time, LH_data,'x','LineWidth',2)
legend("Predicted", "Actual")
title("LH")
xlim([1100, 1153])
xticklabels({})

subplot(4,1,3)
plot(t, E2,'-', time, E2_data,'x','LineWidth',2)
legend("Predicted", "Actual")
title("E_2")
xlim([1100, 1153])
xticklabels({})

subplot(4,1,4)
plot(t, P4, '--',time, P4_data,'x','LineWidth',2)
legend("Predicted", "Actual")
title("P_4")
xlim([1100, 1153])
xticks(1:5)
xticks(1100:10:1153)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})
xlabel("Time (Days)")
sgtitle('Predicted and Actual Hormonal Levels During the Cycle', 'FontSize', 15)




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
sgtitle('Predicted and Actual Hormonal Levels During the Cycle', 'FontSize', 16)


