% @Description: The main algorithm for the oblique plasma wave model 
% features a hybrid distribution that integrates loss-cone product-bi-kappa, 
% kappa-Maxwellian, and bi-Maxwellian plasma distributions.
% @Filename: main_bopbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-1
% @LastEditors: Bai Wei
% @LastEditTime: 2025-01-23

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

deg = 1.0e-4;
theta = deg*pi/180;
%
sp = 0; % sp=0: 'eig()'; sp=1: sparse 'eigs()'; sp>1: Computing Eigenvalues with 'gpuArray' and 'eig()'.
B0 = 1.0e-5;  % background magnetic field in z direction 
nk = 100;
kk0 = linspace(1e-3,2.5,nk);
EPS0 = 1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modules_path = '../../../src';
addpath(modules_path);
par = importdata('./bopbk.in', ' ', 1); % read input parameters
plasmaParams = getPlasmaParameters(B0,par,J_opt); 
if sum(plasmaParams.index_bm)==0
    solver = @(kx,kz,theta) solver_pbk(kx,kz,theta,plasmaParams,sp,EPS0);
elseif sum(plasmaParams.index_pbk)==0
    solver = @(kx,kz,theta) solver_bm(kx,kz,theta,plasmaParams,sp,EPS0);
else
    solver = @(kx,kz,theta) solver_mixed(kx,kz,theta,plasmaParams,sp,EPS0);
end

params_with_unit;

%% Input parameters
S = plasmaParams.S;
Ns = plasmaParams.Ns;
index_pbk = plasmaParams.index_pbk;
index_bm = plasmaParams.index_bm;
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
% kDs = 1./lambdaDs;
% kn = sqrt(kDs(1)^2 + kDs(2)^2);
% wn = sqrt(wps(1)^2 + wps(2)^2);

% Thermal speed of the sth component is vs = (Tsparallel/ms)^0.5;
vs = (kB*Tsz./ms).^0.5;

% gyroradius of the ion is ai = vi/Omegai*(Tiperp/Tiparallel)^0.5;
ai = vs(2)/wcs(2)*sqrt(Tsx(2)/Tsz(2));

%%
tic;
for ik = 1:nk
    k = kk0(ik)*wps(1)/sqrt(c2);
    kx = k*sin(theta); % perpendicular to the magnetic B0
    kz = k*cos(theta); % parallel to the magnetic B0
    [w,eigVec] = solver(kx,kz,theta);
    ww(ik,:) = w;
    www(ik,1,:) = w;
    %
    if(mod(ik,nk/5) == 0)
        fprintf('Progress: %2.2f%%\n',ik/nk*100);
    end
end
runtime = toc;
display_runtime(runtime);
fprintf('\n=====\n');
rmpath(modules_path); % Remove folders from search path
%
kk = kk0*wps(1)/sqrt(c2);
kxx = k*sin(theta); % perpendicular to the magnetic B0
kzz = k*cos(theta); % parallel to the magnetic B0

%%
betasz = 2*mu0*kB.*ns0.*Tsz./B0^2; % beta_para
betasx = 2*mu0*kB.*ns0.*Tsx./B0^2; % beta_perp
vA = B0/sqrt(mu0*sum(ms.*ns0)); % Alfven speed
%
disp(['betaz = ', num2str(betasz)]);
disp(['betax = ', num2str(betasx)]);
disp(['[Alfven speed: vA/c = ', num2str(vA/sqrt(c2),'%.3e')]);
disp(['omega_pe/omega_ce = ', num2str(wps(1)/wcs(1),'%.3e')]);

%% Save data
% Create a new data folder if it doesn't exist
currentPath = pwd;
foldername = 'output';
run('../../../tools/createDateFile(currentPath,foldername)');
fname = ['./',foldername,'/bopbkData.mat'];
save(fname,'deg','theta','S','par','deg','kk0','kk','ai','ww','www','wps','wcs','lambdaDs','rhocs','Ns','runtime');

%% Load data
foldername = 'output';
fname = ['./',foldername,'/bopbkData.mat'];
load(fname);
real_w = real(ww);
imag_w = imag(ww);

%% plot all roots
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',15);
subplot(2,1,1);
plot(kk0,real_w/wps(1),'k.');
xlabel('ck/\omega_{pe}');
ylabel('\omega_r/\omega_{pe}');
grid on;
% xlim([0,0.37]);
% ylim([0 3]);

subplot(2,1,2);
plot(kk0,imag_w/wps(1),'k.','markersize',10);
xlabel('ck/\omega_{pe}');
ylabel('\gamma/\omega_{pe}');
grid on;
% xlim([0,0.37]);
% ylim([-1, 0.25]);

%% plot the selected roots
rootNo = 1:100;
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',15);
subplot(2,1,1);
plot(kk0,real_w(:,rootNo)/abs(wps(1)),'k.','markersize',15);
xlabel('ck/\omega_{pe}');
ylabel('\omega_r/\omega_{pe}');

grid on;
% xlim([0,0.37]);
ylim([0, 1.0e-2]);
title(['$\theta=$',num2str(deg,'%.2f'),'$^\circ$'  ...
       ', $\kappa_{e\parallel}=$',num2str(par.data(1,8)), ...
       ', $\kappa_{e\perp}=$',num2str(par.data(1,9)), ...
       ', $\kappa_{p\parallel}=$',num2str(par.data(2,8)),...
       ', $\kappa_{p\perp}=$',num2str(par.data(2,9))],'Interpreter','latex');
set(gca,'Fontsize',20);

subplot(2,1,2);
plot(kk0,imag_w(:,rootNo)/abs(wps(1)),'b.-','markersize',15);
xlabel('ck/\omega_{pe}');
ylabel('\gamma/\omega_{pe}');
grid on;
