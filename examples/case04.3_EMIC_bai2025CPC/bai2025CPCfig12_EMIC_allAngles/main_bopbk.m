% @Description: The primary algorithm for the oblique plasma wave model 
% exhibits a hybrid distribution that amalgamates product-bi-kappa (PBK), 
% kappa-Maxwellian (KM), Maxwellian-kappa (MK) and bi-Maxwellian (BM) plasma distributions. 
% @Parameter conditions for PBK types: 
%       kappasz>1/2 and kappasx>1 for PBK type I; 
%       kappasz>3/2 and kappasx>2 for PBK type II.
% @Filename: main_bopbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-01
% @LastEditors: Bai Wei
% @LastEditTime: 2026-04-06
% Note: kappasz>1/2 and kappasx>1 for PBK type I; 
%       kappasz>3/2 and kappasx>2 for PBK type II.


% close all;
% clear;
% clc;

% format long

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
B0 = 100.0e-9;  % background magnetic field in z direction 
%
nth = 60;
theta_deg = linspace(1e-4,90-1e-2,nth); % in degree
theta_rad = theta_deg*pi/180; %  in radian
%
nk = 100;
kk0 = linspace(1e-5,0.08,nk);
EPS0 = 1.0e-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
src_root = '../../../src';
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
ms = plasmaParams.ms; % mass
ns0 = plasmaParams.ns0;

%%%
kDs = 1./lambdaDs;
kn = sqrt(kDs(1)^2 + kDs(2)^2);
wn = sqrt(wps(1)^2 + wps(2)^2);

% Thermal speed of the sth component is vs = (Tsparallel/ms)^0.5;
vs = (kB*Tsz./ms).^0.5;

% gyroradius of the ion is ai = vi/Omegai*(Tiperp/Tiparallel)^0.5;
ai = vs(3)/wcs(3)*sqrt(Tsx(3)/Tsz(3));

%%
norm_omega = abs(wcs(1));

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
    if(mod(ith,nth/5) == 0)
        fprintf('Progress: %2.1f%%\n',ith/nth*100);
    end
    %
    for ik = 1:nk
        k = kk0(ik)/rhocs(1);
        kx = k*sin(theta_rad(ith)); % perpendicular to the magnetic B0
        kz = k*cos(theta_rad(ith)); % parallel to the magnetic B0
        % w = double(solver(kx,kz,theta(ith)));
        [w,eigVec] = solver(kx,kz,theta_rad(ith));
        ww(ik,ith,:) = w(:);
    end
end
runtime = toc;
display_runtime(runtime);
fprintf('=====\n');
rmpath(modules_path); % Remove folders from search path
%
kk = kk0/rhocs(1);
% kxx = k*sin(theta); % perpendicular to the magnetic B0
% kzz = k*cos(theta); % parallel to the magnetic B0

%%
betasz = 2*mu0*kB.*ns0.*Tsz./B0^2; % beta_para
betasx = 2*mu0*kB.*ns0.*Tsx./B0^2; % beta_perp
vA = B0./sqrt(mu0*sum(ms.*ns0)); % Alfven speed
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
save(fname,'par','theta_deg','theta_rad','kk0','kk','ww','wps','wcs','lambdaDs','rhocs','Ns','runtime');

%% Load data
foldername = 'output';
fname = ['./',foldername,'/bopbkData.mat'];
load(fname);
real_w = real(ww);
imag_w = imag(ww);
%
[KK0,THETA] = meshgrid(kk0,theta_rad);
KKX0 = KK0.*sin(THETA);
KKZ0= KK0.*cos(THETA);

%% Growth rate and frequency contour plot
ww_contour = ww;
index1 = imag(ww_contour) > 1.0e1; %2e-2;
index2 = imag(ww_contour) < 1.0e-13; %1.0e0; % 0.04;
% index1 = real(ww_contour/wps(1)) > 1.0e-1; %2e-2;
% index2 = real(ww_contour/wps(1)) < 1.0e-1; %1.0e0; % 0.04;
index = index1 | index2;
%
ww_contour(index) = NaN + 1i*NaN;
wr_contour = real(ww_contour);
wi_contour = imag(ww_contour);
% clc
num_start = 1;
num_end = 1;


figure('unit','normalized','Position',[0.01 0.45 0.3 0.65],'DefaultAxesFontSize',20);
subplot(2,1,1)
for ii=num_start:num_end
    contourf(KK0,THETA*90/(pi/2),wr_contour(:,:,ii)','LineColor','none');
    hold on;
end
shading interp
colormap('jet');
colorbar;
xlabel('k\rho_{cp}');
ylabel('\theta^{\circ}');
% xlim([1e-3 1e5]);
% ylim([0 1]);

%
subplot(2,1,2);
for ii=num_start:num_end
    % contourf(KKX0,KKZ0,imag_w(:,:,ii)'/abs(wcs(1)),'LineColor','none');
    contourf(KK0,THETA*90/(pi/2),wi_contour(:,:,ii)',30);
    hold on;
end
shading interp
colormap('jet');
colorbar;
xlabel('k\rho_{cp}');
ylabel('\theta^{\circ}');
% xlim([1e-3 1e5]);
% ylim([0 0.7]);

%% contour plots 
close all;
% h=figure('unit','normalized','Position',[0.01 0.45 0.3 0.65],'DefaultAxesFontSize',20);
figure;
subplot(2,1,1)
for ii=1:10
    contourf(KK0,THETA*180/pi,real_w(:,:,ii)','LineColor','none');
    % surfc(THETA*180/pi,KK0,real_w(:,:,ii)'/abs(wcs(1)));
    hold on;
    %scatter(rkk,rth,[],rww/abs(wcs(1)),'filled','MarkerFaceAlpha',0.6);
end
shading interp
colormap('jet');
colorbar;
xlabel('\theta[^\circ]');
ylabel('k\rho_{cp}');
% xlim([1e-3 1e5]);
% ylim([0 1]);

%
subplot(2,1,2);
for ii=1:1
    contourf(KK0,THETA*360/pi,imag_w(:,:,ii)','LineColor','none');
    % surfc(THETA*180/pi,KK0,imag_w(:,:,ii)'/abs(wcs(1)));
    hold on;
end
shading interp
colormap('jet');
colorbar;
xlabel('\theta[^\circ]');
ylabel('k\rho_{cp}');
% xlim([1e-3 1e5]);
% ylim([0 0.7]);


%% plot the growth rates for the selected angle 
close all;
h=figure('unit','normalized','Position',[0.01 0.45 0.4 0.5],'DefaultAxesFontSize',20);

iTH = 30;
for rootNo = 1:100
subplot(2,1,1);
plot(kk*rhocs(1),real_w(:,iTH,rootNo),'k.','markersize',8);
grid on;
xlabel('k\rho_{cp}');
ylabel('\omega_r/\omega_{cp}');
% xlim([0 0.3]);
ylim([0 1]);
hold on;

subplot(2,1,2);
plot(kk*rhocs(1),imag_w(:,iTH,rootNo),'b.-','markersize',8);
grid on;
xlabel('k\rho_{cp}');
ylabel('\gamma/\omega_{cp}');
% xlim([0 0.3]);
ylim([-0.1,0.1]);
hold on;

% pause(1);
end
