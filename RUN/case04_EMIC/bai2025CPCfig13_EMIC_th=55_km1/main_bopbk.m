% @Description: The main algorithm for the oblique plasma wave model 
% features a hybrid distribution that integrates loss-cone product-bi-kappa, 
% kappa-Maxwellian, and bi-Maxwellian plasma distributions.
% Ref.: BO-KM: Wei Bai et al., A comprehensive solver for dispersion relation 
% of obliquely propagating waves in magnetized multi-species plasma with 
% anisotropic drift kappa-Maxwellian distribution, Computer Physics 
% Communications 307(2025) 109434, https://doi.org/10.1016/j.cpc.2024.109434
% @Filename: main_bopbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-1
% @LastEditors: Bai Wei
% @LastEditTime: 2025.11.16

close all;
clear;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = 8.1;

deg = 55.0;
%
theta = deg*pi/180;
sp = 0; % sp=0: 'eig()'; sp=1: sparse 'eigs()'; sp>1: Computing Eigenvalues with 'gpuArray' and 'eig()'.
B0 = 100.0e-9;  % background magnetic field in z direction 
nk = 100;
kk0 = linspace(1e-5,0.08,nk);

EPS0 = 1e-6; % for integrals S1,S2,...,S9 when \lambda_s \rightarrow 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modules_path = '../../../src';
addpath(modules_path);
par = importdata('./bopbk.in', ' ', 1); % read input parameters
[S,Ns,index_pbk,index_bm,kappasz,kappasx,vtsz,vtsx,Tsz,Tsx,sgms,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaParameters(B0,par); 
if sum(index_bm)==0
    solver = @(kx,kz,theta) solver_pbk(kx,kz,theta,B0,par,sp,EPS0);
elseif sum(index_pbk)==0
    solver = @(kx,kz,theta) solver_bm(J,kx,kz,theta,B0,par,sp,EPS0);
else
    solver = @(kx,kz,theta) solver_mixed(J,kx,kz,theta,B0,par,sp,EPS0);
end

params_with_unit;


%%
tic;
for ik = 1:nk
    k = kk0(ik)/rhocs(1);
    kx = k*sin(theta); % perpendicular to the magnetic B0
    kz = k*cos(theta); % parallel to the magnetic B0
    [w,eigVec] = solver(kx,kz,theta);
    ww(ik,:) = w;
    www(ik,1,:) = w;
    if(mod(ik,nk/5) == 0)
        fprintf('The case completion: %2.1f%%\n',ik/nk*100);
    end
end
rmpath(modules_path); % Remove folders from search path
kk = kk0/rhocs(1);
kxx = k*sin(theta); % perpendicular to the magnetic B0
kzz = k*cos(theta); % parallel to the magnetic B0

runtime = toc;
disp(['Time elapsed: ', num2str(runtime/60), ' minutes']);

%%
for s=1:S
    ms(s)=par.data(s,2)*me; % mass
    ns0(s)=par.data(s,3); % desity unit: m^-3
    Tsz(s) = par.data(s,4)*qe/kB; % parallel temperature, unit: eV -> K
    Tsx(s) = par.data(s,5)*qe/kB; % perp temperature, unit: eV -> K
end

%
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
plot(kk0,real_w(:,:)/abs(wcs(1)),'k.');
ylabel('\omega_r/\omega_{ce}');
xlabel('k\rho_{ce}');
grid on;
ylim([0 1]);

subplot(2,1,2);
plot(kk0,imag_w(:,:)/abs(wcs(1)),'k.','markersize',10);
ylabel('\gamma/\omega_{ce}');
xlabel('k\rho_{ce}');
grid on;
% xlim([0,0.37]);
ylim([-0.05,0.05]);

%% plot the selected roots
rootNo = 1:400;
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',15);
subplot(2,1,1);
plot(kk0,real_w(:,rootNo)/abs(wcs(1)),'k.','markersize',15);
ylabel('\omega_r/\omega_{ce}');
xlabel('k/\omega_{ce}');
% xlim([1e-3 1e5]);
ylim([0 1]);
grid on;
hold on;
title(['$\theta=$',num2str(deg,'%.1f'),'$^\circ$', ...
    ' ($\kappa_{e}=$', num2str(par.data(1,8)),')'], ... 
    'Interpreter','latex','FontSize',20);

subplot(2,1,2);
plot(kk0,imag_w(:,rootNo)/abs(wcs(1)),'b.-','markersize',15);
ylabel('\gamma/\omega_{ce}');
xlabel('k/\omega_{ce}');
grid on;
hold on;
ylim([-0.05,0.05]);

