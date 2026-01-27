% @Description: The main algorithm for the oblique plasma wave model 
% features a hybrid distribution that integrates loss-cone product-bi-kappa, 
% kappa-Maxwellian, and bi-Maxwellian plasma distributions.
% @Filename: main_bokm.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-1
% @LastEditors: Bai Wei
% @LastEditTime: 2025.11.17

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
sp = 0; % sp=0: 'eig()'; sp=1: sparse 'eigs()'; sp>1: Computing Eigenvalues with 'gpuArray' and 'eig()'.
B0 = 100.0e-9;  % background magnetic field in z direction 
%
nth = 60;
deg = linspace(1e-4,90-1e-2,nth); % in degree
theta = deg*pi/180;
%
nk = 100;
kk0 = linspace(1e-5,0.08,nk);
EPS0 = 1.0e-2;
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
tic;
for ith = 1:nth
    for ik = 1:nk
        k = kk0(ik)/rhocs(1);
        kx = k*sin(theta(ith)); % perpendicular to the magnetic B0
        kz = k*cos(theta(ith)); % parallel to the magnetic B0
        % w = double(solver(kx,kz,theta(ith)));
        [w,eigVec] = solver(kx,kz,theta(ith));
        ww(ith,ik,:) = w(:);
    end
    if(mod(ith,nth/5) == 0)
        fprintf('Progress: %2.1f%%\n',ith/nth*100);
    end
end
runtime = toc;
display_runtime(runtime);
fprintf('\n=====\n');
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
save(fname,'par','deg','theta','kk0','kk','ww','wps','wcs','lambdaDs','rhocs','Ns','runtime');

%% Load data
foldername = 'output';
fname = ['./',foldername,'/bopbkData.mat'];
load(fname);
real_w = real(ww);
imag_w = imag(ww);
%
[THETA,KK0] = meshgrid(theta,kk0);
KKX0 = KK0.*sin(THETA);
KKZ0= KK0.*cos(THETA);

%% Growth rate and frequency contour plot
ww_contour = ww;
index1 = imag(ww_contour/wps(1)) > 1.0e1; %2e-2;
index2 = imag(ww_contour/wps(1)) < 1.0e-7; %1.0e0; % 0.04;
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
    contourf(KK0,THETA*90/(pi/2),wr_contour(:,:,ii)'/abs(wcs(1)),'LineColor','none');
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
    contourf(KK0,THETA*90/(pi/2),wi_contour(:,:,ii)'/abs(wcs(1)),30);
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
    contourf(THETA*180/pi,KK0,real_w(:,:,ii)'/abs(wcs(1)),'LineColor','none');
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
    contourf(THETA*360/pi,KK0,imag_w(:,:,ii)'/abs(wcs(1)),'LineColor','none');
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
plot(kk*rhocs(1),real_w(iTH,:,rootNo)/abs(wcs(1)),'k.','markersize',8);
grid on;
xlabel('k\rho_{cp}');
ylabel('\omega_r/\omega_{cp}');
% xlim([0 0.3]);
ylim([0 1]);
hold on;

subplot(2,1,2);
plot(kk*rhocs(1),imag_w(iTH,:,rootNo)/abs(wcs(1)),'b.-','markersize',8);
grid on;
xlabel('k\rho_{cp}');
ylabel('\gamma/\omega_{cp}');
% xlim([0 0.3]);
ylim([-0.05,0.05]);
hold on;

% pause(1);
end

