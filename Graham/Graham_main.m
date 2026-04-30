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
[t,Y] = ode15s(@(t,y) Graham(t,y,graham_params,9.6377), tspan,init); 

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

figure(3);
tcl = tiledlayout(3, 1);

nexttile(tcl)
line1 = plot(t, PH, '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Predicted');
xticklabels({})
xlim([1100, 1153])
title('Follicular Phase');

nexttile(tcl)
line2 = plot(t, OM, '-', 'Color', '#5F2399', 'LineWidth', 2);
xticklabels({})
xlim([1100, 1153])
title('Ovulatory Phase');

nexttile(tcl)
line3 = plot(t, LB, '-', 'Color', '#5F2399', 'LineWidth', 2);
xticklabels({})
xlim([1100, 1153])
title('Luteal Phase');
xticks(1100:10:1153)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70"})

fontsize(gcf, 12, "points")
hL.FontSize = 14;
hL.Layout.Tile = 'east';
sgtitle('Ovulatory Cycle', 'FontSize', 16)


%% Lots of E2

% sequences of parameter values to test 
E0_list = linspace(1,120,10);

tic

for j = 1:length(E0_list)

    % run simulation
    [t,Y] = ode15s(@(t,y) Graham(t,y,graham_params, E0_list(j)), tspan, init);

    % Preallocate after first run
    if j == 1
        Nt = length(t);

        fshp_vals = zeros(length(E0_list), Nt);
        fsh_vals  = zeros(length(E0_list), Nt);
        lhp_vals  = zeros(length(E0_list), Nt);
        lh_vals   = zeros(length(E0_list), Nt);
        ph_vals   = zeros(length(E0_list), Nt);
        om_vals   = zeros(length(E0_list), Nt);
        lb_vals   = zeros(length(E0_list), Nt);
        e2_vals   = zeros(length(E0_list), Nt);
        p4_vals   = zeros(length(E0_list), Nt);
    end

    % Store full trajectories (transpose to make row)
    fshp_vals(j,:) = Y(:,1)';
    fsh_vals(j,:)  = Y(:,2)';
    lhp_vals(j,:)  = Y(:,3)';
    lh_vals(j,:)   = Y(:,4)';
    ph_vals(j,:)   = Y(:,5)';
    om_vals(j,:)   = Y(:,6)';
    lb_vals(j,:)   = Y(:,7)';
    e2_vals(j,:)   = Y(:,9)';
    p4_vals(j,:)   = Y(:,10)';

end

toc


% timing
elapsedTime = toc;  % stop timer
fprintf('Loop runtime: %.2f seconds\n', elapsedTime);
%% plots
% plotting


figure(4);
tcl = tiledlayout(4, 1);

nexttile(tcl)
hold on
line1 = plot(t, fsh_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t-10, fsh_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Normal E_2');
line13 = plot(t, fsh_vals(7,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t, fsh_vals(10,:), ':','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line11 = plot(time, FSH_data, 'x', 'Color', '#3A5ABD','LineWidth', 2, 'DisplayName', 'Data') ;
hold off
xticklabels({})
xlim([1100, 1200])
title('FSH (IU/L)');

nexttile(tcl)
hold on
line2 = plot(t, lh_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line22 = plot(t-10, lh_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Normal E_2');
line23 = plot(t, lh_vals(7,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line24 = plot(t, lh_vals(10,:), ':','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line22 = plot(time, LH_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xticklabels({})
xlim([1100, 1200])
title('LH (IU/L)');

nexttile(tcl)
hold on
line1 = plot(t, e2_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t-10, e2_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Normal E_2');
line13 = plot(t, e2_vals(7,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t, e2_vals(10,:), ':','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line33 = plot(time, E2_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xticklabels({})
xlim([1100, 1200])
title('E_2 (pg/mL)');

nexttile(tcl)
hold on
line1 = plot(t, p4_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t-10, p4_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Normal E_2');
line13 = plot(t, p4_vals(7,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t, p4_vals(10,:), ':','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
%line44 = plot(time, P4_data, 'x', 'Color', '#3A5ABD', 'LineWidth', 2, 'DisplayName','Data');
hold off
xlim([1100, 1200])
xticks(1100:10:1200)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})
xlabel("Time (Days)")
title('P_4 (pg/mL)');

%hL = legend([line1, line12, line13, line14]);
%hL.FontSize = 14;
%hL.Layout.Tile = 'east';
fontsize(gcf, 12, "points")
title(tcl, 'Hormonal Levels During the Cycle', 'FontSize', 16)
subtitle(tcl, 'Varying E_2 Levels', 'FontSize', 12)

figure(5);
tcl = tiledlayout(3, 1);

nexttile(tcl)
hold on
line1 = plot(t, ph_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t-10, ph_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Normal E_2');
line13 = plot(t, ph_vals(7,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t, ph_vals(10,:), ':','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
xticklabels({})
xlim([1100, 1200])
title('Follicular Phase');

nexttile(tcl)
hold on
line1 = plot(t, om_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t-10, om_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Normal E_2');
line13 = plot(t, om_vals(7,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t, om_vals(10,:), ':','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
xticklabels({})
xlim([1100, 1200])
title('Ovulatory Phase');

nexttile(tcl)
hold on
line1 = plot(t, lb_vals(1,:), '-','Color', '#5F2399','LineWidth',2, 'DisplayName','Low E_2');
line12 = plot(t-10, lb_vals(3,:), '--','Color', '#BF40BF','LineWidth',2, 'DisplayName','Normal E_2');
line13 = plot(t, lb_vals(7,:), '-.','Color', '#593163','LineWidth',2, 'DisplayName','High E_2');
line14 = plot(t, lb_vals(10,:), ':','Color', '#480082','LineWidth',2, 'DisplayName','Very High E_2');
xticklabels({})
xlim([1100, 1200])
title('Luteal Phase');
xticks(1100:10:1200)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})

fontsize(gcf, 12, "points")
hL = legend([line1, line12, line13, line14]);
hL.FontSize = 14;
hL.Layout.Tile = 'east';
title(tcl, 'Ovulatory Cycle', 'FontSize', 16)
subtitle(tcl, 'Varying E_2 Levels', 'FontSize', 12)

%% cleaned plotting %%

figure(5);
tcl = tiledlayout(3,1);

% Indices of simulations to plot
idx = [1 3 5 10];

% Labels
labels = {'Low E_2','Normal E_2','High E_2','Very High E_2'};

% Data groups
data_groups = {ph_vals, om_vals, lb_vals};
titles = {'Follicular Phase','Ovulatory Phase','Luteal Phase'};

for k = 1:3
    nexttile
    hold on
    
    for i = 1:length(idx)
        plot(t, data_groups{k}(idx(i),:), ...
            'LineWidth',2, ...
            'DisplayName', labels{i});
    end
    
    xlim([1100 1153])
    title(titles{k})
    
    if k < 3
        xticklabels({})
    else
        xticks(1100:10:1153)
        xticklabels({"0","10","20","30","40","50","60","70"})
    end
end

% Single legend for whole layout
hL = legend;
hL.Layout.Tile = 'east';
hL.FontSize = 14;

sgtitle('Ovulatory Cycle','FontSize',16)
set(gcf,'DefaultAxesFontSize',12)