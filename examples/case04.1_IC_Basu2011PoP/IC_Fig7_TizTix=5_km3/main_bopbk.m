% @Description: Oblique plasma wave solver using a hybrid distribution.
% Includes: loss-cone product-bi-kappa (PBK), kappa-Maxwellian (KM), 
% Maxwellian-kappa (MK), and bi-Maxwellian (BM). 
% @Parameter conditions for PBK types: 
%       kappasz>1/2 and kappasx>1 for PBK type I; 
%       kappasz>3/2 and kappasx>2 for PBK type II.
% Ref. B. Basu and N. J. Grossbard, Ion-cyclotron instability in current-carrying Lorentzian 
% (kappa) and Maxwellian plasmas with anisotropic temperatures: 
% A comparative study, Physics of Plasmas 18, 092106 (2011).
% https://doi.org/10.1063/1.3632974
% @Filename: main_bopbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-01
% @LastEditors: Bai Wei
% @LastEditTime: 2026-01-25

% close all;
% clear;
% clc;

%% input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J=2.1:  Huba2009, 2-pole
% J=2.2:  Martin1979, 2-pole, J=2, I=3
% J=3:    Martin1980, 3-pole, J=3, I=3
% J=4.1:  Martin1980, 4-pole, J=4, I=5 
% J=4.2:  new calculation, J=4, I=5
% J=8.1:  Ronnmark1982, 8-pole for Z function
% J=8.2:  J=8, I=8
% J=8.3:  J=8, I=10
% J=8.4:  optimized J=8 pole from Xie 2024
% J=10.1  J=10 (I=12,K=8) pole from Xie 2024
% J=10.2  J=10 (I=14,K=6) pole from Xie 2024
% J=12.1: J=12; I=16; 2014;
% J=12.2: J=12; I=16; 2018
% J=12.3: J=12; I=12;
% J=16.1: J=16; I=18; 
% J=16.2: J=16 (I=16,K=20) pole from Xie 2024 
% J=16.3: J=16 (J=24,I=8) pole from Xie 2024
% J=20:   J=20 (I=23,K=17) pole from Xie 2024
% J=24.1: (J=24,I=24) pole from BO code
% J=24.2: J=24 (I=24,K=24) pole from Xie 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_opt = 8.1;

%% Use 1 for PBK type I or 2 for PBK type II.
pbk_type = 1; 

%%
sp = 0; % sp=0: 'eig()'; sp=1: sparse 'eigs()'; sp>1: Computing Eigenvalues with 'gpuArray' and 'eig()'.
B0 = 1.0e-3; % background magnetic field in z direction 
%
nk = 100;
nth = 1;
%
kpara0 = 0.08;
kperp0_max = 3.0;
kkperp0 = linspace(1e-3,kperp0_max,nk);
theta = atan(kkperp0/kpara0);
kk0 = kpara0./cos(theta);

EPS0 = 1.0e-2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
src_root = '../../../src';  %
addpath(fullfile(src_root, 'common'));
if pbk_type == 1
    modules_path = fullfile(src_root, 'pbk_type1');  % Type I
    disp("Type I PBK");
elseif pbk_type == 2
    modules_path = fullfile(src_root, 'pbk_type2');  % Type II
    disp("Type II PBK");
else 
    error("Invalid pbk_type value %d. Use 1 for PBK type I or 2 for PBK type II.", pbk_type);
end
addpath(modules_path);
%
par = importdata('./bopbk.in', ' ', 1); % read input parameters
if pbk_type == 1 
    plasmaParams = getPlasmaParameters_pbk1(B0,par,J_opt); 
elseif pbk_type == 2
    plasmaParams = getPlasmaParameters_pbk2(B0,par,J_opt);
end
%
params_with_unit;

%% Input parameters
S = plasmaParams.S;
Ns = plasmaParams.Ns;
pbk_mask = plasmaParams.is_pbk;
mk_mask = plasmaParams.is_mk;
kappasz = plasmaParams.kappasz;
kappasx = plasmaParams.kappasx;
vtsz = plasmaParams.vtsz;
vtsx = plasmaParams.vtsx;
Tsz = plasmaParams.Tsz;
Tsx = plasmaParams.Tsx;
sgms = plasmaParams.sgms;
wps = plasmaParams.wps;
wcs = plasmaParams.wcs;
us0 = plasmaParams.us0;
rhocs = plasmaParams.rhocs;
lambdaDs = plasmaParams.lambdaDs;
ms = plasmaParams.ms;
ns0 = plasmaParams.ns0;

%
me = ms(1);
mi = ms(2);
ne = ns0(1);
ni = ns0(2);

vA_p = B0/sqrt(mu0*mi*ni); % Alfven speed of proton
vA = B0/sqrt(mu0*sum(ms.*ns0)); % Alfven speed

%%
norm_omega = wcs(2);

%% solver
if pbk_type==1
    if sum(mk_mask)==0
        solver = @(kx,kz,theta) solver_pbk1(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0);
    elseif sum(pbk_mask)==0
        solver = @(kx,kz,theta) solver_mk1(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0);
    else
        solver = @(kx,kz,theta) solver_mixed_pbk1(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0);
    end
elseif pbk_type==2
     if sum(mk_mask)==0
        solver = @(kx,kz,theta) solver_pbk2(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0);
    elseif sum(pbk_mask)==0
        solver = @(kx,kz,theta) solver_mk2(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0);
    else
        solver = @(kx,kz,theta) solver_mixed_pbk2(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0);
     end
end

%%
tic;
ith = 1;
for ik = 1:nk
    k = kk0(ik)/(rhocs(2)/sqrt(2));
    kx = k*sin(theta(ik)); % perpendicular to the magnetic B0
    kz = k*cos(theta(ik)); % parallel to the magnetic B0
    [w,eigVec] = solver(kx,kz,theta(ik));
    ww(ik,:) = w;
    www(ik,1,:) = w;
    if(mod(ik,nk/5) == 0)
        fprintf('Progress: %2.1f%%\n',ik/nk*100);
    end
end
runtime = toc;
display_runtime(runtime);
fprintf('=====\n');
rmpath(modules_path); % Remove folders from search path
%
kk = kk0/(rhocs(2)/sqrt(2));
kkx = kk.*sin(theta); % perpendicular to the magnetic B0
kkz = kk.*cos(theta); % parallel to the magnetic B0

%%
betasz = 2*mu0*kB.*ns0.*Tsz./B0^2; % beta_para
betasx = 2*mu0*kB.*ns0.*Tsx./B0^2; % beta_perp
vA = B0/sqrt(mu0*sum(ms.*ns0)); % Alfven speed
%
disp(['betaz = ', num2str(betasz)]);
disp(['betax = ', num2str(betasx)]);
disp(['Alfven speed: vA/c = ', num2str(vA/sqrt(c2))]);

%% Save data
% Create a new data folder if it doesn't exist
currentPath = pwd;
foldername = 'output';
run('../../../tools/createDateFile(currentPath,foldername)');
fname = ['./',foldername,'/bopbkData.mat'];
save(fname,'par','theta','kpara0','kkperp0','kk0',...
    'ww','www','wps','wcs','lambdaDs','rhocs','Ns','runtime');

%% Load data
foldername = 'output';
fname = ['./',foldername,'/bopbkData.mat'];
load(fname);
real_w = real(ww);
imag_w = imag(ww);

%% plot all roots
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',15);
subplot(2,1,1);
plot(kkperp0,real_w,'k.');
xlabel('$k_{\perp}\rho_i/\sqrt{2}$','Interpreter','latex');
ylabel('$\omega_r/\Omega_{i}$','Interpreter','latex');
grid on;
% xlim([0,0.37]);
ylim([0,2]);

subplot(2,1,2);
plot(kkperp0,imag_w,'k.','markersize',10);
xlabel('$k_{\perp}\rho_i/\sqrt{2}$','Interpreter','latex');
ylabel('$\gamma/\Omega_{i}$','Interpreter','latex');
grid on;
% xlim([0,0.37]);
% ylim([0,0.8]);

%% plot the selected roots
rootNo = 1:10;
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',15);
subplot(2,1,1);
plot(kkperp0,real_w(:,rootNo),'k.','markersize',15);
xlabel('$k_{\perp}\rho_i/\sqrt{2}$','Interpreter','latex');
ylabel('$\omega_r/\Omega_{i}$','Interpreter','latex');
grid on;
% xlim([0,0.37]);
% ylim([-100,100]);
title(['$k_{\parallel}\rho_i/\sqrt{2}=$',num2str(kpara0,'%.1f'), ...
    ' ($\kappa_{e}=$', num2str(par.data(1,8)), ...
    ', $\kappa_{p}=$', num2str(par.data(2,8)),')'], ... 
    'Interpreter','latex','FontSize',20);

subplot(2,1,2);
plot(kkperp0,imag_w(:,rootNo),'b.-','markersize',15);
xlabel('$k_{\perp}\rho_i/\sqrt{2}$','Interpreter','latex');
ylabel('$\gamma/\Omega_{i}$','Interpreter','latex');
grid on;
% ylim([0,10]);
