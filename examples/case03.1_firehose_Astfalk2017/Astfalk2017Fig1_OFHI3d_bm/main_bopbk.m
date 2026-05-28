% @Description: The main algorithm for the oblique plasma wave model 
% features a hybrid distribution that integrates product-bi-kappa (PBK), 
% kappa-Maxwellian（KM), Maxwellian-kappa (MK)  and bi-Maxwellian (BM) plasma distributions.
% @Parameter conditions for PBK types: 
%       kappasz>1/2 and kappasx>1 for PBK type I; 
%       kappasz>3/2 and kappasx>2 for PBK type II.
% Ref. Astfalk, P., and F. Jenko (2017), LEOPARD: A grid-based dispersion 
% relation solver for arbitrary gyrotropic distributions, 
% J. Geophys. Res. Space Physics, 122, 89–101, 
% doi:10.1002/2016JA023522.
% @Filename: main_bopbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-01
% @Modification: 
%   Bai Wei, 2026.04.07: Updated BM to MK (retains BM features).

% clear;
% close all;

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
nth = 170;
% nth = 10;
deg = linspace(1e-2,90-1e-2,nth);
theta = deg*pi/180;
sp = 0; % sp=0: eig(); sp=1: sparse eigs(); sp>=2: gpuArray
B0 = 1.0e-1; % background magnetic field in z direction 

nk = 160;
% nk = 12;
kk0 = linspace(1e-3,0.8,nk);
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
for ith = 1:nth
    for ik = 1:nk
        k = kk0(ik)/di;
        kx = k*sin(theta(ith)); % perpendicular to the magnetic B0
        kz = k*cos(theta(ith)); % parallel to the magnetic B0
        w = double(solver(kx,kz,theta(ith)));
        ww(ik,ith,:) = w(:);
    end
    %
    if(mod(ith,nth/5) == 0)
        fprintf('Progress: %2.1f%%\n',ith/nth*100);
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
vA_p = B0./sqrt(mu0*ms(2)*ns0(2)); % Alfven speed of proton
vA = B0./sqrt(mu0*sum(ms.*ns0)); % Alfven speed
%
disp(['betasz = ', num2str(betasz)]);
disp(['betasx = ', num2str(betasx)]);
disp(['[Alfven speed]/c = ', num2str(vA/sqrt(c2),'%.3e')]);

%% Save data
% Create a new data folder if it doesn't exist
currentPath = pwd;
foldername = 'output';
run('../../../tools/createDateFile(currentPath,foldername)');
fname = ['./',foldername,'/Astfalk2017Fig1_OFHI3d_bm.mat'];
save(fname,'par','deg','theta','kk0','kk','ww','wps','wcs','lambdaDs','rhocs','Ns','runtime');

%% Load data
foldername = 'output';
fname = ['./',foldername,'/Astfalk2017Fig1_OFHI3d_bm.mat'];
load(fname);
real_w = real(ww);
imag_w = imag(ww);
[THETA,KK0] = meshgrid(theta,kk0);

%% Growth rate and frequency contour plot
ww_contour = ww;
index1 = imag(ww_contour/wcs(2)) < 5.0e-3;
index2 = imag(ww_contour/wcs(2)) > 1.0e1;
index = index1 | index2;
ww_contour(index) = NaN + 1i*NaN;
realw_contour = real(ww_contour);
imagw_contour = imag(ww_contour);
% clc
num_contour = 1;

h=figure('unit','normalized','Position',[0.01 0.2 0.4 0.5],'DefaultAxesFontSize',20);
for ii=1:num_contour
    % contourf(KK0,THETA*180/pi,realw_contour(:,:,ii)/wcs(3));
    contourf(KK0(:,:),THETA(:,:)*180/pi,realw_contour(:,:,ii)/wcs(2),20);
    hold on;
end
shading interp
colormap('jet');
colorbar;
xlabel('kd_i');
ylabel('\theta');
% xlim([1e-3 1e5]);
% ylim([0 1e1]);
title(['$\omega_r/\Omega_{p}$: ', '$\kappa_{e\parallel}=$',num2str(par.data(1,8)), ...
      '$, \kappa_{p\parallel}=$',num2str(par.data(2,8))],'Interpreter','latex');

h=figure('unit','normalized','Position',[0.01 0.2 0.4 0.5],'DefaultAxesFontSize',20);
for ii=1:num_contour
    % contourf(KK0,THETA*180/pi,imagw_contour(:,:,ii)/wcs(3));
    contourf(KK0(:,:),THETA(:,:)*180/pi,imagw_contour(:,:,ii)/wcs(2),20);
    hold on;
end
shading interp
colormap('jet');
colorbar;
xlabel('kd_i');
ylabel('\theta');
% xlim([1e-3 1e5]);
% ylim([0 0.7]);
title(['$\gamma/\Omega_{p}$: ', '$\kappa_{e\parallel}=$',num2str(par.data(1,8)), ...
      '$, \kappa_{p\parallel}=$',num2str(par.data(2,8))],'Interpreter','latex');
