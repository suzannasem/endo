%%%% Incorporating E2 %%%%
%%% Simulation file %%%

%% Initial conditions %%

m00 = 9.0 * 10^5 ;
m10 = 100 ;
m20 = 4.5 * 10^4 ;
k00 = 5 * 10^5 ; 
ka0 = 2 * 10^4 ; 
e00 = 0 ; ef0 = 0; ea0 = 0;
Ra0 = 0; Rb0 = 0; 

init = [m00, m10, m20, k00, ka0, e00, ef0, ea0, Ra0, Rb0] ;

opts = odeset('RelTol',1e-3,'AbsTol',1e-6);

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

% periodic cubic spline interpolation
E2_pp = csape(t, E2_data, 'periodic');
P4_pp = csape(t, P4_data, 'periodic');

% define period of 31 days
T = t(end) - t(1);

% using mod arithmetic to repeat spline 
E2 = @(tq) ppval(E2_pp, mod(tq - t(1), T) + t(1));
P4 = @(tq) ppval(P4_pp, mod(tq - t(1), T) + t(1));

tq = linspace(t(1), 1200, 1200);  % evaluate up to 1200 timesteps

% fitted, periodic splines
E2_vals = E2(tq); 
P4_vals = P4(tq);

%% ODE Simulations %%

tspan = 0:0.5:1200 ;
[t,Y] = ode15s(@(t,y) Variable_E2_fn(t,y, 1e-6, 1e-5, E2, P4, 3e3, 4e3),tspan,init, opts); % disease free

% different values of K1, K2
[t1,Y1] = ode15s(@(t,y) Variable_E2_fn(t,y, 1e-6, 1e-5, E2, P4, 4e3, 3e3),tspan,init, opts); % slight malfunction
[t2,Y2] = ode15s(@(t,y) Variable_E2_fn(t,y, 1e-6, 1e-5, E2, P4, 5e3, 3e3),tspan,init, opts); % early disease
[t3,Y3] = ode15s(@(t,y) Variable_E2_fn(t,y, 1e-6, 1e-5, E2, P4, 6e3, 3e3),tspan,init, opts); % disease

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

% disease free profiles

figure(101)

subplot(2,2,1)
plot(t, k0, '--', t, ka, '-', 'LineWidth', 2)
legend("K_0", "K_A")
ylabel("Cells / mL")
xlim([1010, 1072])
%ylim([0,400e3])
xticks([])
title("Natural Killer Cells")

subplot(2,2,2)
plot(t, e0, 'LineWidth', 2)
ylabel("Cells / mL")
xlim([1010, 1072])
xticks([])
%ylim([0, 15e3])
title("Eutopic Endometrial Cells")

subplot(2,2,3)
semilogy(t, m0, '--', t, m1, '-', t, m2, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
xlabel("Time in days")
ylabel("Cells / mL (log scale)")
xlim([1010, 1072])
xticks(1010:10:1072)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})
xlabel("Time (Days)")
title("Macrophage")

subplot(2,2,4)
plot(t, ea, '-', t, ef, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1010, 1072])
xticks(1010:10:1072)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})
xlabel("Time (Days)")
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells")

% Figure 1

figure(102)

subplot(2,2,1)
plot(t3, k03, '--', t3, ka3, '-', 'LineWidth', 2)
legend("K_0", "K_A")
ylabel("Cells / mL")
xlim([1010, 1072])
xticks([])
%ylim([0,400e3])
title("Natural Killer Cells")

subplot(2,2,2)
plot(t3, e03, 'LineWidth', 2)
ylabel("Cells / mL")
xlim([1010, 1072])
xticks([])
%ylim([0, 15e3])
title("Eutopic Endometrial Cells")

subplot(2,2,3)
semilogy(t3, m03, '--', t3, m13, '-', t3, m23, '-.', 'LineWidth', 2)
legend("M_0", "M_1", "M_2")
ylabel("Cells / mL (log scale)")
xlim([1010, 1072])
xticks(1010:10:1072)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})
xlabel("Time (Days)")
title("Macrophage")

subplot(2,2,4)
plot(t3, ea3, '-', t3, ef3, '-', 'LineWidth', 2)
legend("E_A", "E_F")
xticks(1010:10:1072)
xticklabels({"0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"})
xlabel("Time (Days)")
ylabel("Cells / mL")
xlim([1010, 1072])
%ylim([0,0.7e3])
title("Ectopic Endometrial Cells")

%figure(202)

%plot(t, Ra, '-', t, Rb, '-', 'LineWidth', 2)
%legend("ER_\alpha", "ER_\beta")
%xlabel("Time in days")
%ylabel("Receptors / mL")
%xlim([1010, 1072])
%ylim([0,0.7e3])
%title("Estrogen Receptors on Macrophages")

%subplot(2,1,2)
%plot(t, Da, '-', t, Dab, '-', 'LineWidth', 2)
%legend("D_a", "D_{ab}")
%xlabel("Time in days")
%ylabel("Dimers / mL")
%xlim([1010, 1072])
%ylim([0,0.7e3])
%title("Estrogen Receptor Dimers on Macrophages")

% ERb/ERa ratio (important in distinguishing disease states)

%ratio = Rb ./ Ra ;

%figure(203)
%plot(t, ratio, '-', 'LineWidth',2)
%xlabel('Time in Days')
%title('ER_\beta to ER_\alpha Ratio')
%xlim([1010, 1072])
%ylim([0,8e3])

figure(637)

plot(t, Rb, '--', 'LineWidth', 2)
hold on
plot(t1, Rb1, '-', 'LineWidth', 2)
plot(t2, Rb2, '-.', 'LineWidth', 2)
plot(t3, Rb3, ':', 'LineWidth', 2)
hold off
legend("Disease-Free", "Borderline", "Early Disease", "Late Disease")
xlabel("Time in days")
ylabel("Receptors / mL")
xlim([1010, 1072])
ylim([0,1])
title("ER_\beta Proportion")

figure(638)

subplot(2,1,1)
plot(t, ea, '--', 'LineWidth', 2)
hold on
plot(t1, ea1, '-', 'LineWidth', 2)
hold off
legend("Disease-Free", "Borderline")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1010, 1072])
%ylim([0,0.04])
title("Lesion Cell Count")

subplot(2,1,2)
plot(t2, ea2, '-.', 'LineWidth', 2)
hold on
plot(t3, ea3, ':', 'LineWidth', 2)
hold off
legend("Early Disease", "Late Disease")
xlabel("Time in days")
ylabel("Cells / mL")
xlim([1010, 1072])
ylim([1e7, 1e9])
