% Gradient_plots
% Author: Yunchen Xiao
% This MATLAB file generates the plots of all the 
% spatial/temporal gradients at different 
% locations of the 1D spatial domain (boundaries excluded), averaged across
% all 200 datasets at different CVs.

% Environment setting
clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 18)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

x11 = linspace(0,1,80);
x11_res = x11(2:79);

% Read in all the data
grad_true = readtable('Mean reference data grads.txt');

tc_temporal = readtable('Mean temporal grads tc.txt');
dn_spatial = readtable('Mean spatial grads dn.txt');
ga_spatial = readtable('Mean spatial grads gamma.txt');
rn_spatial = readtable('Mean spatial grads rn.txt');

ecm_temporal = readtable('Mean temporal grads ecm.txt');
ita_spatial = readtable('Mean spatial grads ita.txt');

mde_temporal = readtable('Mean temporal grads mde.txt');
dm_spatial = readtable('Mean spatial grads dm.txt');
alpha_spatial = readtable('Mean spatial grads alpha.txt');

tc_temp_true = table2array(grad_true(1,2:79));
dn_spat_true = table2array(grad_true(2,2:79));
ga_spat_true = table2array(grad_true(3,2:79));
rn_spat_true = table2array(grad_true(4,2:79));
ecm_temp_true = table2array(grad_true(5,2:79));
ita_spat_true = table2array(grad_true(6,2:79));
mde_temp_true = table2array(grad_true(7,2:79));
dm_spat_true = table2array(grad_true(8,2:79));
alpha_spat_true = table2array(grad_true(9,2:79));

tc_temp_cv1 = table2array(tc_temporal(1,2:79));
tc_temp_cv2 = table2array(tc_temporal(2,2:79));
tc_temp_cv3 = table2array(tc_temporal(3,2:79));
tc_temp_cv4 = table2array(tc_temporal(4,2:79));
tc_temp_cv5 = table2array(tc_temporal(5,2:79));

dn_spat_cv1 = table2array(dn_spatial(1,2:79));
dn_spat_cv2 = table2array(dn_spatial(2,2:79));
dn_spat_cv3 = table2array(dn_spatial(3,2:79));
dn_spat_cv4 = table2array(dn_spatial(4,2:79));
dn_spat_cv5 = table2array(dn_spatial(5,2:79));

ga_spat_cv1 = table2array(ga_spatial(1,2:79));
ga_spat_cv2 = table2array(ga_spatial(2,2:79));
ga_spat_cv3 = table2array(ga_spatial(3,2:79));
ga_spat_cv4 = table2array(ga_spatial(4,2:79));
ga_spat_cv5 = table2array(ga_spatial(5,2:79));

rn_spat_cv1 = table2array(rn_spatial(1,2:79));
rn_spat_cv2 = table2array(rn_spatial(2,2:79));
rn_spat_cv3 = table2array(rn_spatial(3,2:79));
rn_spat_cv4 = table2array(rn_spatial(4,2:79));
rn_spat_cv5 = table2array(rn_spatial(5,2:79));

ecm_temp_cv1 = table2array(ecm_temporal(1,2:79));
ecm_temp_cv2 = table2array(ecm_temporal(2,2:79));
ecm_temp_cv3 = table2array(ecm_temporal(3,2:79));
ecm_temp_cv4 = table2array(ecm_temporal(4,2:79));
ecm_temp_cv5 = table2array(ecm_temporal(5,2:79));

ita_spat_cv1 = table2array(ita_spatial(1,2:79));
ita_spat_cv2 = table2array(ita_spatial(2,2:79));
ita_spat_cv3 = table2array(ita_spatial(3,2:79));
ita_spat_cv4 = table2array(ita_spatial(4,2:79));
ita_spat_cv5 = table2array(ita_spatial(5,2:79));

mde_temp_cv1 = table2array(mde_temporal(1,2:79));
mde_temp_cv2 = table2array(mde_temporal(2,2:79));
mde_temp_cv3 = table2array(mde_temporal(3,2:79));
mde_temp_cv4 = table2array(mde_temporal(4,2:79));
mde_temp_cv5 = table2array(mde_temporal(5,2:79));

dm_spat_cv1 = table2array(dm_spatial(1,2:79));
dm_spat_cv2 = table2array(dm_spatial(2,2:79));
dm_spat_cv3 = table2array(dm_spatial(3,2:79));
dm_spat_cv4 = table2array(dm_spatial(4,2:79));
dm_spat_cv5 = table2array(dm_spatial(5,2:79));

alpha_spat_cv1 = table2array(alpha_spatial(1,2:79));
alpha_spat_cv2 = table2array(alpha_spatial(2,2:79));
alpha_spat_cv3 = table2array(alpha_spatial(3,2:79));
alpha_spat_cv4 = table2array(alpha_spatial(4,2:79));
alpha_spat_cv5 = table2array(alpha_spatial(5,2:79));

% Temporal gradients of tumour cell density (dn/dt)
figure
plot(x11_res, tc_temp_true, 'r--')
hold on;
plot(x11_res, tc_temp_cv1,'k-')
plot(x11_res, tc_temp_cv2,'g-')
plot(x11_res, tc_temp_cv1,'b-')
plot(x11_res, tc_temp_cv1,'c-')
plot(x11_res, tc_temp_cv1,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{\frac{\partial n}{\partial t}}$')

% Diffusion term of tumour cells
figure
plot(x11_res, dn_spat_true, 'r--')
hold on;
plot(x11_res, dn_spat_cv1,'k-')
plot(x11_res, dn_spat_cv2,'g-')
plot(x11_res, dn_spat_cv3,'b-')
plot(x11_res, dn_spat_cv4,'c-')
plot(x11_res, dn_spat_cv5,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{\frac{\partial^{2}n}{\partial x^{2}}}$')

% Haptotaxis term
figure
plot(x11_res, ga_spat_true, 'r--')
hold on;
plot(x11_res, ga_spat_cv1,'k-')
plot(x11_res, ga_spat_cv2,'g-')
plot(x11_res, ga_spat_cv3,'b-')
plot(x11_res, ga_spat_cv4,'c-')
plot(x11_res, ga_spat_cv5,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{\frac{\partial}{\partial x}\left(n\frac{\partial f}{\partial x}\right)}$')

% Logistic growth of tumour cells
figure
plot(x11_res, rn_spat_true, 'r--')
hold on;
plot(x11_res, rn_spat_cv1,'k-')
plot(x11_res, rn_spat_cv2,'g-')
plot(x11_res, rn_spat_cv3,'b-')
plot(x11_res, rn_spat_cv4,'c-')
plot(x11_res, rn_spat_cv5,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{n\left(1-n-f\right)}$')

% Temporal gradients of ECM density (df\dt)
figure
plot(x11_res, ecm_temp_true, 'r--')
hold on;
plot(x11_res, ecm_temp_cv1,'k-')
plot(x11_res, ecm_temp_cv2,'g-')
plot(x11_res, ecm_temp_cv3,'b-')
plot(x11_res, ecm_temp_cv4,'c-')
plot(x11_res, ecm_temp_cv5,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{\frac{\partial f}{\partial t}}$')

% Decay of ECM.
figure
plot(x11_res, ita_spat_true, 'r--')
hold on;
plot(x11_res, ita_spat_cv1,'k-')
plot(x11_res, ita_spat_cv2,'g-')
plot(x11_res, ita_spat_cv3,'b-')
plot(x11_res, ita_spat_cv4,'c-')
plot(x11_res, ita_spat_cv5,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{mf}$')

% Temporal gradients of MDE. (dm/dt)
figure
plot(x11_res, mde_temp_true, 'r--')
hold on;
plot(x11_res, mde_temp_cv1,'k-')
plot(x11_res, mde_temp_cv2,'g-')
plot(x11_res, mde_temp_cv3,'b-')
plot(x11_res, mde_temp_cv4,'c-')
plot(x11_res, mde_temp_cv5,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{\frac{\partial m}{\partial t}}$')

% Diffusion term of MDE.
figure
plot(x11_res, dm_spat_true, 'r--')
hold on;
plot(x11_res, dm_spat_cv1,'k-')
plot(x11_res, dm_spat_cv2,'g-')
plot(x11_res, dm_spat_cv3,'b-')
plot(x11_res, dm_spat_cv4,'c-')
plot(x11_res, dm_spat_cv5,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{\frac{\partial^{2}m}{\partial x^{2}}}$')

% Growth/proliferation of MDE.
figure
plot(x11_res, alpha_spat_true, 'r--')
hold on;
plot(x11_res, alpha_spat_cv1,'k-')
plot(x11_res, alpha_spat_cv2,'g-')
plot(x11_res, alpha_spat_cv3,'b-')
plot(x11_res, alpha_spat_cv4,'c-')
plot(x11_res, alpha_spat_cv5,'m-')
%xlim([0 0.11])
%ylim([0 0.004])
xlabel('Spatial domain')
ylabel('Mean estimated gradients')
title('$\hat{n}$')
