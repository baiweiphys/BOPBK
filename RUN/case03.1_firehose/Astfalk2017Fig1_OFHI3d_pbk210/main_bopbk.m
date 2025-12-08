% @Description: The main algorithm for the oblique plasma wave model 
% features a hybrid distribution that integrates loss-cone product-bi-kappa, 
% kappa-Maxwellian, and bi-Maxwellian plasma distributions.
% Ref. Astfalk, P., and F. Jenko (2017), LEOPARD: A grid-based dispersion 
% relation solver for arbitrary gyrotropic distributions, 
% J. Geophys. Res. Space Physics, 122, 89–101, 
% doi:10.1002/2016JA023522.
% @Filename: main_bopbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-1
% @LastEditors: Bai Wei
% @LastEditTime: 2025.11.28

clear;
close all;

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

nth = 170;
deg = linspace(1e-2,90-1e-2,nth);
theta = deg*pi/180;
sp = 1; % sp=0: eig(); sp=1: sparse eigs(); sp>=2: gpuArray
B0 = 1.0e-1; % background magnetic field in z direction 

nk = 160;
kk0 = linspace(1e-3,0.8,nk);
EPS0 = 1.0e-2; 

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
di = vA_p/wcs(2); % ion inertial length di=vA/Omega_i

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
        fprintf('The case completion: %2.1f%%\n',ith/nth*100);
    end
end

rmpath(modules_path); % Remove folders from search path
kk = kk0/di;
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
betasz = 2*mu0*kB.*ns0'.*Tsz./B0^2; % beta_para
betasx = 2*mu0*kB.*ns0'.*Tsx./B0^2; % beta_perp
vA_p = B0./sqrt(mu0*ms(2)*ns0(2)); % Alfven speed of proton
vA = B0./sqrt(mu0*sum(ms'.*ns0)); % Alfven speed
%
disp(['betasz = ', num2str(betasz)]);
disp(['betasx = ', num2str(betasx)]);
disp(['[Alfven speed]/c = ', num2str(vA/sqrt(c2),'%.3e')]);

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

[THETA,KK0] = meshgrid(theta,kk0);
% KKX = KK0.*sin(THETA);
% KKZ= KK0.*cos(THETA);

%% Growth rate and frequency contour plot
ww_contour = ww;
index1 = imag(ww_contour/wcs(2)) < 1.0e-3;
index2 = imag(ww_contour/wcs(2)) > 2.0e1;
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
title(['$\omega_r/\Omega_{p}$: ', ...
      '$\kappa_{e\parallel}=$',num2str(par.data(1,8)), ...
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
title(['$\gamma/\Omega_{p}$: ', ...
      '$\kappa_{e\parallel}=$',num2str(par.data(1,8)), ...
      '$, \kappa_{p\parallel}=$',num2str(par.data(2,8))],'Interpreter','latex');
