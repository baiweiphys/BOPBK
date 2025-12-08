% @Description: The main algorithm for the oblique plasma wave model 
% features a hybrid distribution that integrates loss-cone product-bi-kappa, 
% kappa-Maxwellian, and bi-Maxwellian plasma distributions.
% Ref. B. Basu and N. J. Grossbard, Ion-cyclotron instability in current-carrying Lorentzian 
% (kappa) and Maxwellian plasmas with anisotropic temperatures: 
% A comparative study, Physics of Plasmas 18, 092106 (2011).
% https://doi.org/10.1063/1.3632974
% @Filename: main_bopbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-1
% @LastEditors: Bai Wei
% @LastEditTime: 2025.02.19

close all;
clear;

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
sp = 0; % sp=0: 'eig()'; sp=1: sparse 'eigs()'; sp>1: Computing Eigenvalues with 'gpuArray' and 'eig()'.
B0 = 1.0e-3; % background magnetic field in z direction 
%
nk = 150;
%
kpara0 = 0.08;
kperp0_max = 2.0;
kkperp0 = linspace(1e-3,kperp0_max,nk);
theta = atan(kkperp0/kpara0);
kk0 = kpara0./cos(theta);
EPS0 = 1e-2; 

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
ms = par.data(:,2)*me; % mass
ns0 = par.data(:,3); % desity unit: m^-3
me = ms(1);
mi = ms(2);
ne = ns0(1);
ni = ns0(2);

vA_p = B0/sqrt(mu0*mi*ni); % Alfven speed of proton
vA = B0/sqrt(mu0*sum(ms.*ns0)); % Alfven speed

%%
tic;
for ik = 1:nk
    k = kk0(ik)/(rhocs(2)/sqrt(2));
    kx = k*sin(theta(ik)); % perpendicular to the magnetic B0
    kz = k*cos(theta(ik)); % parallel to the magnetic B0
    [w,eigVec] = solver(kx,kz,theta(ik));
    ww(ik,:) = w;
    www(ik,1,:) = w;
    if(mod(ik,nk/5) == 0)
        fprintf('The case completion: %2.1f%%\n',ik/nk*100);
    end
end
rmpath(modules_path); % Remove folders from search path
kk = kk0/(rhocs(2)/sqrt(2));
kkx = kk.*sin(theta); % perpendicular to the magnetic B0
kkz = kk.*cos(theta); % parallel to the magnetic B0

runtime = toc;
disp(['Time elapsed: ', num2str(runtime/60), ' minutes']);

%%
for s=1:S
    ms(s)=par.data(s,2)*me; % mass
    ns0(s)=par.data(s,3); % desity unit: m^-3
    Tsz(s) = par.data(s,4)*qe/kB; % parallel temperature, unit: eV -> K
    Tsx(s) = par.data(s,5)*qe/kB; % perp temperature, unit: eV -> K
end

betasz = 2*mu0*kB.*ns0'.*Tsz./B0^2; % beta_para
betasx = 2*mu0*kB.*ns0'.*Tsx./B0^2; % beta_perp
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
plot(kkperp0,real_w/wcs(2),'k.');
xlabel('$k_{\perp}\rho_i/\sqrt{2}$','Interpreter','latex');
ylabel('$\omega_r/\Omega_{i}$','Interpreter','latex');
grid on;
% xlim([0,0.37]);
ylim([0,2]);

subplot(2,1,2);
plot(kkperp0,imag_w/wcs(2),'k.','markersize',10);
xlabel('$k_{\perp}\rho_i/\sqrt{2}$','Interpreter','latex');
ylabel('$\gamma/\Omega_{i}$','Interpreter','latex');
grid on;
% xlim([0,0.37]);
% ylim([0,0.8]);

%% plot the selected roots
rootNo = 1:10;
h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.6],'DefaultAxesFontSize',15);
subplot(2,1,1);
plot(kkperp0,real_w(:,rootNo)/wcs(2),'k.','markersize',15);
xlabel('$k_{\perp}\rho_i/\sqrt{2}$','Interpreter','latex');
ylabel('$\omega_r/\Omega_{i}$','Interpreter','latex');
grid on;
% xlim([0,0.37]);
% ylim([-100,100]);
title(['$k_{\parallel}\rho_i/\sqrt{2}=$',num2str(kpara0,'%.2f'), ...
    ' ($\kappa_{e}=$', num2str(par.data(1,8)), ...
    ', $\kappa_{p}=$', num2str(par.data(2,8)),')'], ... 
    'Interpreter','latex','FontSize',20);

subplot(2,1,2);
plot(kkperp0,imag_w(:,rootNo)/wcs(2),'b.-','markersize',15);
xlabel('$k_{\perp}\rho_i/\sqrt{2}$','Interpreter','latex');
ylabel('$\gamma/\Omega_{i}$','Interpreter','latex');
grid on;
% ylim([0,10]);
