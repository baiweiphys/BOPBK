% @Description: The main algorithm for the oblique plasma wave model 
% features a hybrid distribution that integrates product-bi-kappa (PBK), 
% kappa-Maxwellian（KM), Maxwellian-kappa (MK)  and bi-Maxwellian (BM) plasma distributions.
% @Parameter conditions for PBK types: 
%       kappasz>1/2 and kappasx>1 for PBK type I; 
%       kappasz>3/2 and kappasx>2 for PBK type II.
% @Filename: main_bopbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-01
% @Modification: 
%   Bai Wei, 2026.04.07: Updated BM to MK (retains BM features).


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
deg = 45.0;
theta = deg*pi/180;
sp = 0; % sp=0: 'eig()'; sp=1: sparse 'eigs()'; sp>1: Computing Eigenvalues with 'gpuArray' and 'eig()'.
B0 = 0.1;   % background magnetic field in z direction 
nk = 100;
kk0 = linspace(1e-3,0.8,nk);
% kk0 = logspace(-4,-1,nk);
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
di = vA_p/wcs(2); % ion inertial length di=vA/Omega_i

%%
norm_omega = 1;

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
for ik = 1:nk
    k = kk0(ik)/di;
    kx = k*sin(theta); % perpendicular to the magnetic B0
    kz = k*cos(theta); % parallel to the magnetic B0
    % [eigenVector,w] = double(solver(kx,kz,theta));
    [w,eigVec] = solver(kx,kz,theta);
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
kk = kk0/di;
kxx = k*sin(theta); % perpendicular to the magnetic B0
kzz = k*cos(theta); % parallel to the magnetic B0

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
save(fname,'par','deg','theta','kk0','kk','ww','www','wps','wcs','lambdaDs','rhocs','Ns','runtime');

%% Load data
foldername = 'output';
fname = ['./',foldername,'/bopbkData.mat'];
load(fname);
real_w = real(ww);
imag_w = imag(ww);

%% plot all roots
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',15);
subplot(2,1,1);
plot(kk0,real_w/wcs(2),'k.');
xlabel('kd_i');
ylabel('\omega_r/\Omega_{i}');
grid on;
% xlim([0,0.37]);
% ylim([0 3]);

subplot(2,1,2);
plot(kk0,imag_w/wcs(2),'k.','markersize',10);
xlabel('kd_i');
ylabel('\gamma/\Omega_{i}');
grid on;
% xlim([0,0.37]);
% ylim([-1, 0.25]);

%% plot the selected roots
rootNo = 1:100;
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',15);
subplot(2,1,1);
plot(kk0,real_w(:,rootNo)/wcs(2),'k.','markersize',15);
ylabel('\omega_r/\Omega_{i}');
xlabel('kd_i');
grid on;
% xlim([0,0.37]);
% ylim([0 1.4]);
title(['$\theta=$',num2str(deg,'%.1f'),'$^\circ$', ...
    ' ($\kappa_{e}=$', num2str(par.data(1,8)), ...
    ', $\kappa_{p}=$', num2str(par.data(2,8)),')'], ... 
    'Interpreter','latex','FontSize',20);

subplot(2,1,2);
plot(kk0,imag_w(:,rootNo)/wcs(2),'b.-','markersize',15);
xlabel('kd_i');
ylabel('\gamma/\Omega_{i}');
grid on;
ylim([-1e-1, 1e-1]);
