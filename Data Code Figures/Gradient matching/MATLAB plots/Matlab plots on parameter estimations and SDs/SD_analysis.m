% SD_analysis.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of SD (standard deviation) estimations 
% under different CVs (Coefficient of variations) obtained using three different 
% method: true standard deviations (standard deviations of the parameter
% estimates obtained in 200 perturbed datasets.); analytical standard
% deviations (standard deviations obtained using inverse hessian);
% bootstrap standard deviations.

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

% Read in all the data
std_true = readtable('Sample standard deviations.txt');
std_hess = readtable('Hessian standard deviations.txt');
std_boots = readtable('Bootstrap standard deviations.txt');

% SDs of dn
dn_std_samp = table2array(std_true(1:5,2));
dn_std_hess = table2array(std_hess(1:5,2));
dn_std_boots = table2array(std_boots(1:5,2));

% SDs of gamma
ga_std_samp = table2array(std_true(1:5,3));
ga_std_hess = table2array(std_hess(1:5,3));
ga_std_boots = table2array(std_boots(1:5,3));

% SDs of rn
rn_std_samp = table2array(std_true(1:5,4));
rn_std_hess = table2array(std_hess(1:5,4));
rn_std_boots = table2array(std_boots(1:5,4));

% SDs of eta
eta_std_samp = table2array(std_true(1:5,5));
eta_std_hess = table2array(std_hess(1:5,5));
eta_std_boots = table2array(std_boots(1:5,5));

% SDs of dm
dm_std_samp = table2array(std_true(1:5,6));
dm_std_hess = table2array(std_hess(1:5,6));
dm_std_boots = table2array(std_boots(1:5,6));

% SDs of alpha
alpha_std_samp = table2array(std_true(1:5,7));
alpha_std_hess = table2array(std_hess(1:5,7));
alpha_std_boots = table2array(std_boots(1:5,7));

x = [0.01 0.025 0.05 0.075 0.1];

% Plot of SDs of dn obtained using the three different methods mentioned
% above
figure
plot(x, dn_std_samp, 'r--')
hold on;
plot(x, dn_std_hess,'k-')
plot(x, dn_std_boots,'b-')
plot(0.01, dn_std_samp(1), 'xr', 'markersize', 20)
plot(0.025, dn_std_samp(2), 'xr', 'markersize', 20)
plot(0.05, dn_std_samp(3), 'xr', 'markersize', 20)
plot(0.075, dn_std_samp(4), 'xr', 'markersize', 20)
plot(0.1, dn_std_samp(5), 'xr', 'markersize', 20)
plot(0.01, dn_std_hess(1), 'xk', 'markersize', 20)
plot(0.025, dn_std_hess(2), 'xk', 'markersize', 20)
plot(0.05, dn_std_hess(3), 'xk', 'markersize', 20)
plot(0.075, dn_std_hess(4), 'xk', 'markersize', 20)
plot(0.1, dn_std_hess(5), 'xk', 'markersize', 20)
plot(0.01, dn_std_boots(1), 'xb', 'markersize', 20)
plot(0.025, dn_std_boots(2), 'xb', 'markersize', 20)
plot(0.05, dn_std_boots(3), 'xb', 'markersize', 20)
plot(0.075, dn_std_boots(4), 'xb', 'markersize', 20)
plot(0.1, dn_std_boots(5), 'xb', 'markersize', 20)
xlim([0 0.11])
ylim([0 0.004])
xlabel('Level of perturbations')
ylabel('Standard deviations of $d_{n}$')
title('$d_{n}$')

% Plot of SDs of gamma obtained using the three different methods mentioned
% above
figure
plot(x, ga_std_samp, 'r--')
hold on;
plot(x, ga_std_hess,'k-')
plot(x, ga_std_boots,'b-')
plot(0.01, ga_std_samp(1), 'xr', 'markersize', 20)
plot(0.025, ga_std_samp(2), 'xr', 'markersize', 20)
plot(0.05, ga_std_samp(3), 'xr', 'markersize', 20)
plot(0.075, ga_std_samp(4), 'xr', 'markersize', 20)
plot(0.1, ga_std_samp(5), 'xr', 'markersize', 20)
plot(0.01, ga_std_hess(1), 'xk', 'markersize', 20)
plot(0.025, ga_std_hess(2), 'xk', 'markersize', 20)
plot(0.05, ga_std_hess(3), 'xk', 'markersize', 20)
plot(0.075, ga_std_hess(4), 'xk', 'markersize', 20)
plot(0.1, ga_std_hess(5), 'xk', 'markersize', 20)
plot(0.01, ga_std_boots(1), 'xb', 'markersize', 20)
plot(0.025, ga_std_boots(2), 'xb', 'markersize', 20)
plot(0.05, ga_std_boots(3), 'xb', 'markersize', 20)
plot(0.075, ga_std_boots(4), 'xb', 'markersize', 20)
plot(0.1, ga_std_boots(5), 'xb', 'markersize', 20)
xlim([0 0.11])
ylim([0 0.007])
xlabel('Level of perturbations')
ylabel('Standard deviations of $\gamma$')
title('$\gamma$')

% Plot of SDs of rn obtained using the three different methods mentioned
% above
figure
plot(x, rn_std_samp, 'r--')
hold on;
plot(x, rn_std_hess,'k-')
plot(x, rn_std_boots,'b-')
plot(0.01, rn_std_samp(1), 'xr', 'markersize', 20)
plot(0.025, rn_std_samp(2), 'xr', 'markersize', 20)
plot(0.05, rn_std_samp(3), 'xr', 'markersize', 20)
plot(0.075, rn_std_samp(4), 'xr', 'markersize', 20)
plot(0.1, rn_std_samp(5), 'xr', 'markersize', 20)
plot(0.01, rn_std_hess(1), 'xk', 'markersize', 20)
plot(0.025, rn_std_hess(2), 'xk', 'markersize', 20)
plot(0.05, rn_std_hess(3), 'xk', 'markersize', 20)
plot(0.075, rn_std_hess(4), 'xk', 'markersize', 20)
plot(0.1, rn_std_hess(5), 'xk', 'markersize', 20)
plot(0.01, rn_std_boots(1), 'xb', 'markersize', 20)
plot(0.025, rn_std_boots(2), 'xb', 'markersize', 20)
plot(0.05, rn_std_boots(3), 'xb', 'markersize', 20)
plot(0.075, rn_std_boots(4), 'xb', 'markersize', 20)
plot(0.1, rn_std_boots(5), 'xb', 'markersize', 20)
xlim([0 0.11])
ylim([0 0.8])
xlabel('Level of perturbations')
ylabel('Standard deviations of $r_{n}$')
title('$r_{n}$')

% Plot of SDs of eta obtained using the three different methods mentioned
% above
figure
plot(x, eta_std_samp, 'r--')
hold on;
plot(x, eta_std_hess,'k-')
plot(x, eta_std_boots,'b-')
plot(0.01, eta_std_samp(1), 'xr', 'markersize', 20)
plot(0.025, eta_std_samp(2), 'xr', 'markersize', 20)
plot(0.05, eta_std_samp(3), 'xr', 'markersize', 20)
plot(0.075, eta_std_samp(4), 'xr', 'markersize', 20)
plot(0.1, eta_std_samp(5), 'xr', 'markersize', 20)
plot(0.01, eta_std_hess(1), 'xk', 'markersize', 20)
plot(0.025, eta_std_hess(2), 'xk', 'markersize', 20)
plot(0.05, eta_std_hess(3), 'xk', 'markersize', 20)
plot(0.075, eta_std_hess(4), 'xk', 'markersize', 20)
plot(0.1, eta_std_hess(5), 'xk', 'markersize', 20)
plot(0.01, eta_std_boots(1), 'xb', 'markersize', 20)
plot(0.025, eta_std_boots(2), 'xb', 'markersize', 20)
plot(0.05, eta_std_boots(3), 'xb', 'markersize', 20)
plot(0.075, eta_std_boots(4), 'xb', 'markersize', 20)
plot(0.1, eta_std_boots(5), 'xb', 'markersize', 20)
xlim([0 0.11])
ylim([0 2.5])
xlabel('Level of perturbations')
ylabel('Standard deviations of $\eta$')
title('$\eta$')

% Plot of SDs of dm obtained using the three different methods mentioned
% above
figure
plot(x, dm_std_samp, 'r--')
hold on;
plot(x, dm_std_hess,'k-')
plot(x, dm_std_boots,'b-')
plot(0.01, dm_std_samp(1), 'xr', 'markersize', 20)
plot(0.025, dm_std_samp(2), 'xr', 'markersize', 20)
plot(0.05, dm_std_samp(3), 'xr', 'markersize', 20)
plot(0.075, dm_std_samp(4), 'xr', 'markersize', 20)
plot(0.1, dm_std_samp(5), 'xr', 'markersize', 20)
plot(0.01, dm_std_hess(1), 'xk', 'markersize', 20)
plot(0.025, dm_std_hess(2), 'xk', 'markersize', 20)
plot(0.05, dm_std_hess(3), 'xk', 'markersize', 20)
plot(0.075, dm_std_hess(4), 'xk', 'markersize', 20)
plot(0.1, dm_std_hess(5), 'xk', 'markersize', 20)
plot(0.01, dm_std_boots(1), 'xb', 'markersize', 20)
plot(0.025, dm_std_boots(2), 'xb', 'markersize', 20)
plot(0.05, dm_std_boots(3), 'xb', 'markersize', 20)
plot(0.075, dm_std_boots(4), 'xb', 'markersize', 20)
plot(0.1, dm_std_boots(5), 'xb', 'markersize', 20)
xlim([0 0.11])
ylim([0 0.01])
xlabel('Level of perturbations')
ylabel('Standard deviations of $d_{m}$')
title('$d_{m}$')

% Plot of SDs of alpha obtained using the three different methods mentioned
% above
figure
plot(x, alpha_std_samp, 'r--')
hold on;
plot(x, alpha_std_hess,'k-')
plot(x, alpha_std_boots,'b-')
plot(0.01, alpha_std_samp(1), 'xr', 'markersize', 20)
plot(0.025, alpha_std_samp(2), 'xr', 'markersize', 20)
plot(0.05, alpha_std_samp(3), 'xr', 'markersize', 20)
plot(0.075, alpha_std_samp(4), 'xr', 'markersize', 20)
plot(0.1, alpha_std_samp(5), 'xr', 'markersize', 20)
plot(0.01, alpha_std_hess(1), 'xk', 'markersize', 20)
plot(0.025, alpha_std_hess(2), 'xk', 'markersize', 20)
plot(0.05, alpha_std_hess(3), 'xk', 'markersize', 20)
plot(0.075, alpha_std_hess(4), 'xk', 'markersize', 20)
plot(0.1, alpha_std_hess(5), 'xk', 'markersize', 20)
plot(0.01, alpha_std_boots(1), 'xb', 'markersize', 20)
plot(0.025, alpha_std_boots(2), 'xb', 'markersize', 20)
plot(0.05, alpha_std_boots(3), 'xb', 'markersize', 20)
plot(0.075, alpha_std_boots(4), 'xb', 'markersize', 20)
plot(0.1, alpha_std_boots(5), 'xb', 'markersize', 20)
xlim([0 0.11])
%ylim([0 0.02])
xlabel('Level of perturbations')
ylabel('Standard deviations of $\alpha$')
title('$\alpha$')

