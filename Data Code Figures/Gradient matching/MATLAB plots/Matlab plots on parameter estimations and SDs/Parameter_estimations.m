% Parameter_estimations.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of parameter estimations under
% different CVs (Coefficient of variations), the plots of results from the 
% three sensitivity tests are also presented here.

% Environment settings
clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 18)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

% Read in all the data.
paras_all = readtable('Parameters means at all CV levels.txt');
paras_cons_dndm = readtable('Parameters means at all CV levels dn dm fixed.txt');
paras_cons_dngamma = readtable('Parameters means at all CV levels dn gamma fixed.txt');
paras_cons_dngammarn = readtable('Parameters means at all CV levels dn gamma rn fixed.txt');

% dn estimations
dn_est_all = table2array(paras_all(:,2));

% gamma estimations
gamma_est_all = table2array(paras_all(:,3));
gamma_est_consdndm = table2array(paras_cons_dndm(:,3));

% rn estimations
rn_est_all = table2array(paras_all(:,4));
rn_est_consdndm = table2array(paras_cons_dndm(:,4));
rn_est_consdngamma = table2array(paras_cons_dngamma(:,4));

% eta estimations
eta_est_all = table2array(paras_all(:,5));
eta_est_consdndm = table2array(paras_cons_dndm(:,5));
eta_est_consdngamma = table2array(paras_cons_dngamma(:,5));
eta_est_consdngammarn = table2array(paras_cons_dngammarn(:,5));

% dm estimations
dm_est_all = table2array(paras_all(:,6));
dm_est_consdngamma = table2array(paras_cons_dngamma(:,6));
dm_est_consdngammarn = table2array(paras_cons_dngammarn(:,6));

% alpha estimations
alpha_est_all = table2array(paras_all(:,7));
alpha_est_consdndm = table2array(paras_cons_dndm(:,7));
alpha_est_consdngamma = table2array(paras_cons_dngamma(:,7));
alpha_est_consdngammarn = table2array(paras_cons_dngammarn(:,7));


x = [0.01 0.025 0.05 0.075 0.1];                               

% Plot of dn estimations under different CVs, no parameter being fixed.
figure
yyaxis left
plot(x, dn_est_all, 'k-')
title('$\hat{d_n}$')
hold on;
plot(0.01, dn_est_all(1), 'xk', 'markersize', 20)
plot(0.025, dn_est_all(2), 'xk', 'markersize', 20)
plot(0.05, dn_est_all(3), 'xk', 'markersize', 20)
plot(0.075, dn_est_all(4), 'xk', 'markersize', 20)
plot(0.1, dn_est_all(5), 'xk', 'markersize', 20)
yline(0.01,'r--','Linewidth',3.5);
xlim([0 0.11])
xlabel('Measurement errors CV')
ylim([0.0035 0.011])
ylabel('Mean ($\hat{d_n}$)')
yyaxis right
h = plot(x,dn_est_all,'k-');
delete(h)
ylim([-65 10])
ylabel('Percentage error')

% Plot of gamma estimations under different CVs in different circumstances:
% No parameters being fixed; dn & dm being fixed
figure
yyaxis left
plot(x, gamma_est_all,'k-')
hold on;
plot(x, gamma_est_consdndm, 'b-')
plot(0.01, gamma_est_all(1), 'xk', 'markersize', 20)
plot(0.025, gamma_est_all(2), 'xk', 'markersize', 20)
plot(0.05, gamma_est_all(3), 'xk', 'markersize', 20)
plot(0.075, gamma_est_all(4), 'xk', 'markersize', 20)
plot(0.1, gamma_est_all(5), 'xk', 'markersize', 20)
plot(0.01, gamma_est_consdndm(1), 'xb', 'markersize', 20)
plot(0.025, gamma_est_consdndm(2), 'xb', 'markersize', 20)
plot(0.05, gamma_est_consdndm(3), 'xb', 'markersize', 20)
plot(0.075, gamma_est_consdndm(4), 'xb', 'markersize', 20)
plot(0.1, gamma_est_consdndm(5), 'xb', 'markersize', 20)
xlim([0 0.11])
ylim([0.025 0.06])
yline(0.05, 'r--','Linewidth',3.5);
xlabel('Measurement errors CV')
ylabel('Mean ($\hat{\gamma}$)')
title('$\hat{\gamma}$')
yyaxis right
h = plot(x,dn_est_all,'k-');
delete(h)
ylim([-50 10])
ylabel('Percentage error')

% Plot of rn estimations under different CVs in different circumstances:
% No parameters being fixed; dn & dm being fixed; dn & gamma being fixed
figure
yyaxis left;
plot(x, rn_est_all,'k-')
hold on;
plot(x, rn_est_consdndm, 'b-')
plot(x, rn_est_consdngamma, 'r-')
plot(0.01, rn_est_all(1), 'xk', 'markersize', 20)
plot(0.025, rn_est_all(2), 'xk', 'markersize', 20)
plot(0.05, rn_est_all(3), 'xk', 'markersize', 20)
plot(0.075, rn_est_all(4), 'xk', 'markersize', 20)
plot(0.1, rn_est_all(5), 'xk', 'markersize', 20)
plot(0.01, rn_est_consdndm(1), 'xb', 'markersize', 20)
plot(0.025, rn_est_consdndm(2), 'xb', 'markersize', 20)
plot(0.05, rn_est_consdndm(3), 'xb', 'markersize', 20)
plot(0.075, rn_est_consdndm(4), 'xb', 'markersize', 20)
plot(0.1, rn_est_consdndm(5), 'xb', 'markersize', 20)
plot(0.01, rn_est_consdngamma(1), 'xr', 'markersize', 20)
plot(0.025, rn_est_consdngamma(2), 'xr', 'markersize', 20)
plot(0.05, rn_est_consdngamma(3), 'xr', 'markersize', 20)
plot(0.075, rn_est_consdngamma(4), 'xr', 'markersize', 20)
plot(0.1, rn_est_consdngamma(5), 'xr', 'markersize', 20)
xlim([0 0.11])
ylim([2 16])
yline(5, 'r--','Linewidth',3.5);
xlabel('Measurement errors CV')
ylabel('Mean ($\hat{r_n})$')
title('$\hat{r_n}$')
yyaxis right
h = plot(x,dn_est_all,'k-');
delete(h)
ylim([-60 220])
ylabel('Percentage error')

% Plot of eta estimations under different CVs in different circumstances:
% No parameters being fixed; dn & dm being fixed; dn & gamma being fixed;
% dn & gamma & rn being fixed
figure
yyaxis left;
plot(x, eta_est_all,'k-')
hold on;
plot(x, eta_est_consdndm, 'b-')
plot(x, eta_est_consdngamma, 'r-')
plot(x, eta_est_consdngammarn, 'g-')
plot(0.01, eta_est_all(1), 'xk', 'markersize', 20)
plot(0.025, eta_est_all(2), 'xk', 'markersize', 20)
plot(0.05, eta_est_all(3), 'xk', 'markersize', 20)
plot(0.075, eta_est_all(4), 'xk', 'markersize', 20)
plot(0.1, eta_est_all(5), 'xk', 'markersize', 20)
plot(0.01, eta_est_consdndm(1), 'xb', 'markersize', 20)
plot(0.025, eta_est_consdndm(2), 'xb', 'markersize', 20)
plot(0.05, eta_est_consdndm(3), 'xb', 'markersize', 20)
plot(0.075, eta_est_consdndm(4), 'xb', 'markersize', 20)
plot(0.1, eta_est_consdndm(5), 'xb', 'markersize', 20)
plot(0.01, eta_est_consdngamma(1), 'xr', 'markersize', 20)
plot(0.025, eta_est_consdngamma(2), 'xr', 'markersize', 20)
plot(0.05, eta_est_consdngamma(3), 'xr', 'markersize', 20)
plot(0.075, eta_est_consdngamma(4), 'xr', 'markersize', 20)
plot(0.1, eta_est_consdngamma(5), 'xr', 'markersize', 20)
plot(0.01, eta_est_consdngammarn(1), 'xg', 'markersize', 20)
plot(0.025, eta_est_consdngammarn(2), 'xg', 'markersize', 20)
plot(0.05, eta_est_consdngammarn(3), 'xg', 'markersize', 20)
plot(0.075, eta_est_consdngammarn(4), 'xg', 'markersize', 20)
plot(0.1, eta_est_consdngammarn(5), 'xg', 'markersize', 20)
xlim([0 0.11])
ylim([3 13])
yline(10, 'r--', 'Linewidth', 3.5);
xlabel('Measurement errors CV')
ylabel('Mean ($\hat{\eta}$)')
title('$\hat{\eta}$')
yyaxis right
h = plot(x,dn_est_all,'k-');
delete(h)
ylim([-70 30])
ylabel('Percentage error')

% Plot of dm estimations under different CVs in different circumstances:
% No parameters being fixed; dn & gamma being fixed; dn & gamma & rn being fixed
figure
yyaxis left
plot(x, dm_est_all,'k-')
hold on;
plot(x, dm_est_consdngamma, 'r-')
plot(x, dm_est_consdngammarn, 'g-')
plot(0.01, dm_est_all(1), 'xk', 'markersize', 20)
plot(0.025, dm_est_all(2), 'xk', 'markersize', 20)
plot(0.05, dm_est_all(3), 'xk', 'markersize', 20)
plot(0.075, dm_est_all(4), 'xk', 'markersize', 20)
plot(0.1, dm_est_all(5), 'xk', 'markersize', 20)
plot(0.01, dm_est_consdngamma(1), 'xr', 'markersize', 20)
plot(0.025, dm_est_consdngamma(2), 'xr', 'markersize', 20)
plot(0.05, dm_est_consdngamma(3), 'xr', 'markersize', 20)
plot(0.075, dm_est_consdngamma(4), 'xr', 'markersize', 20)
plot(0.1, dm_est_consdngamma(5), 'xr', 'markersize', 20)
plot(0.01, dm_est_consdngammarn(1), 'xg', 'markersize', 20)
plot(0.025, dm_est_consdngammarn(2), 'xg', 'markersize', 20)
plot(0.05, dm_est_consdngammarn(3), 'xg', 'markersize', 20)
plot(0.075, dm_est_consdngammarn(4), 'xg', 'markersize', 20)
plot(0.1, dm_est_consdngammarn(5), 'xg', 'markersize', 20)
xlim([0 0.11])
ylim([0.0075 0.015]);
yline(0.01, 'r--', 'Linewidth', 3.5);
xlabel('Measurement errors CV')
ylabel('Mean ($\hat{d_m}$)')
title('$\hat{d_m}$')
yyaxis right
h = plot(x,dn_est_all,'k-');
delete(h)
ylim([-25 50])
ylabel('Percentage error')

% Plot of alpha estimations under different CVs in different circumstances:
% No parameters being fixed; dn & dm being fixed; dn & gamma being fixed;
% dn & gamma & rn being fixed
figure
yyaxis left;
plot(x, alpha_est_all,'k-')
hold on;
plot(x, alpha_est_consdndm, 'b-')
plot(x, alpha_est_consdngamma, 'r-')
plot(x, alpha_est_consdngammarn, 'g-')
plot(0.01, alpha_est_all(1), 'xk', 'markersize', 20)
plot(0.025, alpha_est_all(2), 'xk', 'markersize', 20)
plot(0.05, alpha_est_all(3), 'xk', 'markersize', 20)
plot(0.075, alpha_est_all(4), 'xk', 'markersize', 20)
plot(0.1, alpha_est_all(5), 'xk', 'markersize', 20)
plot(0.01, alpha_est_consdndm(1), 'xb', 'markersize', 20)
plot(0.025, alpha_est_consdndm(2), 'xb', 'markersize', 20)
plot(0.05, alpha_est_consdndm(3), 'xb', 'markersize', 20)
plot(0.075, alpha_est_consdndm(4), 'xb', 'markersize', 20)
plot(0.1, alpha_est_consdndm(5), 'xb', 'markersize', 20)
plot(0.01, alpha_est_consdngamma(1), 'xr', 'markersize', 20)
plot(0.025, alpha_est_consdngamma(2), 'xr', 'markersize', 20)
plot(0.05, alpha_est_consdngamma(3), 'xr', 'markersize', 20)
plot(0.075, alpha_est_consdngamma(4), 'xr', 'markersize', 20)
plot(0.1, alpha_est_consdngamma(5), 'xr', 'markersize', 20)
plot(0.01, alpha_est_consdngammarn(1), 'xg', 'markersize', 20)
plot(0.025, alpha_est_consdngammarn(2), 'xg', 'markersize', 20)
plot(0.05, alpha_est_consdngammarn(3), 'xg', 'markersize', 20)
plot(0.075, alpha_est_consdngammarn(4), 'xg', 'markersize', 20)
plot(0.1, alpha_est_consdngammarn(5), 'xg', 'markersize', 20)
xlim([0 0.11])
ylim([0.095 0.12])
yline(0.1,'r--','Linewidth',3.5);
xlabel('Measurement errors CV')
ylabel('Mean ($\hat{\alpha}$)')
title('$\hat{\alpha}$')
yyaxis right
h = plot(x,dn_est_all,'k-');
delete(h)
ylim([-5 20])
ylabel('Percentage error')
