% @Description: Parameter setup script for the benchmark case.

clear all;
clc;

% add path
modules_path = '../../../src';
addpath(modules_path);

params_with_unit;

%%
% read input parameters
par = importdata('./bopbk.in', ' ', 1); % read input parameters
B0 = 9.7251e-5; % B0=5nT for solar wind
J_opt = 8.1;
modules_path = '../../../src';
addpath(modules_path);
plasmaParams = getPlasmaParameters(B0,par,J_opt); 

%% Input parameters
S = plasmaParams.S;
Ns = plasmaParams.Ns;
pbk_mask = plasmaParams.is_pbk;
bm_mask = plasmaParams.is_bm;
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
ns0=par.data(:,3); % desity unit: m^-3
ne = ns0(1);
ni = ns0(2);

%%
for s=1:S
    ms(s)=par.data(s,2)*me; % mass
    ns0(s)=par.data(s,3); % desity unit: m^-3
    Tsz(s) = par.data(s,4)*qe/kB; % parallel temperature, unit: eV -> K
    Tsx(s) = par.data(s,5)*qe/kB; % perp temperature, unit: eV -> K
    vts_paral(s) = sqrt(2*kB*Tsz(s)/ms(s)); 
end

betas_para = 2*mu0*kB.*ns0'.*Tsz./B0^2; % beta_para
betas_perp = 2*mu0*kB.*ns0'.*Tsx./B0^2; % beta_perp

vA_p = B0/sqrt(mu0*ms(2)*ns0(2)); % Alfven speed of proton
vA = B0/sqrt(mu0*sum(ms.*ns0')); % Alfven speed
vte = vts_paral(1)/sqrt(c2); % thermal velocity for electron
vti = vts_paral(2)/sqrt(c2); % thermal velocity for ion

%%
% betai = vtiz^2/vA^2
betai = 2.0;
%
mi = ms(2);
% Tempurature (Unit: eV)
Ti_para = betai*mi*vA_p^2/2.0/kB;

%%
disp(['ne = ', num2str(ne,'%.6e'), ' m^-3']);
disp(['ni = ', num2str(ni,'%.6e'), ' m^-3']);
disp(['beta_para = ', num2str(betas_para)]);
disp(['beta_perp = ', num2str(betas_perp)]);
disp(['Alfven speed of proton: vA_p/c = ', num2str(vA_p/sqrt(c2),'%.4e')]);
disp(['Alfven speed: vA/c = ', num2str(vA/sqrt(c2),'%.4e')]);
disp(['vte/c = ', num2str(vte,'%.4e')]); % thermal velocity for electron
disp(['vti/c = ', num2str(vti,'%.4e')]); % thermal velocity for ion
