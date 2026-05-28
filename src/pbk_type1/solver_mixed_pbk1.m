function [w,eigVec] = solver_mixed_pbk1(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0)
% SOLVER_MIXED To solve for the wave mode roots in the mixed solver. 
% @Description: To solve for the wave mode roots of the oblique plasma wave 
% model using the mixed distribution in both loss-cone bi-Maxwellian and PBK 
% plasmas, the transverse wave number kx and longitudinal wave number kz 
% are specified.
% @Filename: solver_mixed_pbk1.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-8-16
% @Modification: 
%   Bai Wei, 2026.04.01: add norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM structure).
%                        Uses BM properties only if (kappaxs > kappaxs_th); 
%                        otherwise uses MK.
%  Bai Wei, 2026.04.05: Update solver to PBK type 1.


params_with_unit;

% Input parameters
% S = plasmaParams.S;
Ns = plasmaParams.Ns;
J_opt = plasmaParams.J_opt;
pbk_mask = plasmaParams.is_pbk;
mk_mask = plasmaParams.is_mk;
kappasz = plasmaParams.kappasz;
kappasx = plasmaParams.kappasx;
kappasx_th = plasmaParams.kappasx_threshold;
vtsz = plasmaParams.vtsz;
vtsx = plasmaParams.vtsx;
Tsz = plasmaParams.Tsz;
Tsx = plasmaParams.Tsx;
sgms = plasmaParams.sgms;
wps = plasmaParams.wps;
wcs = plasmaParams.wcs;
us0 = plasmaParams.us0;
rhocs = plasmaParams.rhocs;
% is_pbk_vts = plasmaParams.is_pbk_vts;
% lambdaDs = plasmaParams.lambdaDs;

Jpole = floor(J_opt);  % Use the floor function to get the integer part

%
lambdas = 0.5*kx^2*rhocs.^2; % for argument of the modified bessel function

% for PBK distribution 
S_pbk = sum(pbk_mask);
Ns_pbk = Ns(pbk_mask);
kappasz_pbk = kappasz(pbk_mask);
kappasx_pbk = kappasx(pbk_mask);
vtsz_pbk = vtsz(pbk_mask);
vtsx_pbk = vtsx(pbk_mask);
% Tsz_pbk = Tsz(pbk_mask);
% Tsx_pbk = Tsx(pbk_mask);
wps_pbk = wps(pbk_mask);
wcs_pbk = wcs(pbk_mask);
us0_pbk = us0(pbk_mask);
% rhocs_pbk = rhocs(pbk_mask);
% lambdaDs_pbk = lambdaDs(pbk_mask);
lambdas_pbk = lambdas(pbk_mask);
sgms_pbk = sgms(pbk_mask); % 2025.01.22


% for Maxwellian-kappa (MK) plasma
S_mk = sum(mk_mask);
Ns_mk = Ns(mk_mask);
% kappasx_mk = kappasx(index_mk);
vtsz_mk = vtsz(mk_mask);
% vtsx_mk = vtsx(index_mk);
Tsz_mk = Tsz(mk_mask);
Tsx_mk = Tsx(mk_mask);
wps_mk = wps(mk_mask);
wcs_mk = wcs(mk_mask);
us0_mk = us0(mk_mask);
% rhocs_mk = rhocs(index_mk);
% lambdaDs_mk = lambdaDs(index_mk);
lambdas_mk = lambdas(mk_mask);
sgms_mk = sgms(mk_mask); % 2025.01.22

%% for PBK distribution plasma

% common/PBK type 1 or 2
csn = @(s,n) pbk_csn(s,n,kz,kappasz_pbk,vtsz_pbk,wcs_pbk,us0_pbk,norm_omega);

% PBK type 1
b1snl = @(s,n,l) pbk1_b1snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,norm_omega,EPS0);
b2snl = @(s,n,l) pbk1_b2snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,norm_omega,EPS0);
b3snl = @(s,n,l) pbk1_b3snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,norm_omega,EPS0);
b4snl = @(s,n,l) pbk1_b4snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,norm_omega,EPS0);
b5snl = @(s,n,l) pbk1_b5snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,norm_omega,EPS0);
b6snl = @(s,n,l) pbk1_b6snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,norm_omega,EPS0);

% common/PBK type 1 or 2
bx10_pbk = pbk_bx10(wps_pbk,kappasx_pbk,sgms_pbk,norm_omega);
bx11snl = @(s,n,l) pbk_bx11snl(s,n,l,b1snl,wps_pbk,norm_omega);
bx12snl = @(s,n,l) pbk_bx12snl(s,n,l,b2snl,wps_pbk,norm_omega);
bx21snl = @(s,n,l) pbk_bx21snl(s,n,l,b3snl,wps_pbk,norm_omega);
bx22snl = @(s,n,l) pbk_bx22snl(s,n,l,b4snl,wps_pbk,norm_omega);
bx31snl = @(s,n,l) pbk_bx31snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk,norm_omega);
bx32snl = @(s,n,l) pbk_bx32snl(s,n,l,theta,csn,b2snl,wps_pbk,wcs_pbk,norm_omega);
bx33snl = @(s,n,l) pbk_bx33snl(s,n,l,theta,b1snl,wps_pbk,wcs_pbk,norm_omega);

% common/PBK type 1 or 2
by11snl = @(s,n,l) -1*bx21snl(s,n,l);
by12snl = @(s,n,l) -1*bx22snl(s,n,l);
by20_pbk = 1i*epsilon0*sum(wps_pbk.^2)/norm_omega^2; % fixed bug by baiwei, 2025.10.22
by21snl = @(s,n,l) pbk_by21snl(s,n,l,b5snl,wps_pbk,norm_omega);
by22snl = @(s,n,l) pbk_by22snl(s,n,l,b6snl,wps_pbk,norm_omega);
by31snl = @(s,n,l) pbk_by31snl(s,n,l,theta,csn,b3snl,b4snl,wps_pbk,wcs_pbk,norm_omega);
by32snl = @(s,n,l) pbk_by32snl(s,n,l,theta,csn,b4snl,wps_pbk,wcs_pbk,norm_omega);
by33snl = @(s,n,l) pbk_by33snl(s,n,l,theta,b3snl,wps_pbk,wcs_pbk,norm_omega);

% common/PBK type 1 or 2
bz11snl = @(s,n,l) bx31snl(s,n,l);
bz12snl = @(s,n,l) bx32snl(s,n,l);
bz13snl = @(s,n,l) bx33snl(s,n,l);
bz21snl = @(s,n,l) -1*by31snl(s,n,l);
bz22snl = @(s,n,l) -1*by32snl(s,n,l);
bz23snl = @(s,n,l) -1*by33snl(s,n,l);
bz31snl = @(s,n,l) pbk_bz31snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk,norm_omega);
bz32snl = @(s,n,l) pbk_bz32snl(s,n,l,theta,csn,b2snl,wps_pbk,wcs_pbk,norm_omega);
bz33snl = @(s,n,l) pbk_bz33snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk,norm_omega);
bz34snl = @(s,n,l) pbk_bz34snl(s,n,l,theta,b1snl,wps_pbk,wcs_pbk,norm_omega);

%% for Maxwellian-kappa (MK) distribution plasma
% common/PBK type 1 or 2
csnj = @(s,n,jj) mk_csnj(s,n,jj,kz,J_opt,vtsz_mk,wcs_mk,us0_mk,norm_omega);

% PBK type 1
b12snj = @(s,n,jj) mk1_b12snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz_mk,Tsz_mk,Tsx_mk,sgms_mk,wcs_mk,lambdas_mk,norm_omega,EPS0);
b34snj = @(s,n,jj) mk1_b34snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz_mk,Tsz_mk,Tsx_mk,sgms_mk,wcs_mk,lambdas_mk,norm_omega,EPS0);
b56snj = @(s,n,jj) mk1_b56snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz_mk,Tsz_mk,Tsx_mk,sgms_mk,wcs_mk,lambdas_mk,norm_omega,EPS0);

% common/PBK type 1 or 2
bx1_mk = mk_bx1(S_mk,Ns_mk,Jpole,b12snj,csnj,wps_mk,norm_omega);
bx1snj = @(s,n,jj) mk_bx1snj(s,n,jj,b12snj,csnj,wps_mk,norm_omega);
bx2_mk = mk_bx2(S_mk,Ns_mk,Jpole,b34snj,csnj,wps_mk,norm_omega);
bx2snj = @(s,n,jj) mk_bx2snj(s,n,jj,b34snj,csnj,wps_mk,norm_omega);
bx3_mk = mk_bx3(S_mk,Ns_mk,Jpole,theta,b12snj,csnj,wps_mk,norm_omega);
bx3snj = @(s,n,jj) mk_bx3snj(s,n,jj,theta,b12snj,csnj,wps_mk,wcs_mk,norm_omega);

% common/PBK type 1 or 2
by1_mk = mk_by1(S_mk,Ns_mk,Jpole,b34snj,csnj,wps_mk,norm_omega);
by1snj = @(s,n,jj) mk_by1snj(s,n,jj,b34snj,csnj,wps_mk,norm_omega);
by2_mk = mk_by2(S_mk,Ns_mk,Jpole,b56snj,csnj,wps_mk,norm_omega);
by2snj = @(s,n,jj) mk_by2snj(s,n,jj,b56snj,csnj,wps_mk,norm_omega);
by3_mk = mk_by3(S_mk,Ns_mk,Jpole,theta,b34snj,csnj,wps_mk,norm_omega);
by3snj = @(s,n,jj) mk_by3snj(s,n,jj,theta,b34snj,csnj,wps_mk,wcs_mk,norm_omega);

% common/PBK type 1 or 2
bz1_mk = mk_bz1(S_mk,Ns_mk,Jpole,theta,b12snj,csnj,wps_mk,norm_omega);
bz1snj = @(s,n,jj) mk_bz1snj(s,n,jj,theta,b12snj,csnj,wps_mk,wcs_mk,norm_omega);
bz2_mk = mk_bz2(S_mk,Ns_mk,Jpole,theta,b34snj,csnj,wps_mk,norm_omega);
bz2snj = @(s,n,jj) mk_bz2snj(s,n,jj,theta,b34snj,csnj,wps_mk,wcs_mk,norm_omega);
bz3_mk = mk_bz3(S_mk,Ns_mk,Jpole,theta,b12snj,csnj,wps_mk,norm_omega);
bz3snj = @(s,n,jj) mk_bz3snj(s,n,jj,theta,b12snj,csnj,wps_mk,wcs_mk,norm_omega);

%% Step 1: PBK plasma matrix
% PBK type 1
Mx_pbk = Mxy_pbk1_mixed(S_pbk,Ns_pbk,Ns_mk,kappasz_pbk,Jpole,csn, ...
                        bx11snl,bx21snl,bx31snl,bx12snl,bx22snl,bx32snl, ...
                        bx33snl,bx10_pbk,1,5,4,3,5); % Matrix No.1

My_pbk = Mxy_pbk1_mixed(S_pbk,Ns_pbk,Ns_mk,kappasz_pbk,Jpole,csn, ...
                        by11snl,by21snl,by31snl,by12snl,by22snl,by32snl, ...
                        by33snl,by20_pbk,2,5,4,3,4); % Matrix No.2

Mz_pbk = Mz_pbk1_mixed(S_pbk,Ns_pbk,Ns_mk,kappasz_pbk,Jpole,csn, ...
                       bz11snl,bz21snl,bz31snl,bz13snl,bz23snl,bz33snl, ...
                       bz12snl,bz22snl,bz32snl,bz34snl,by20_pbk,3,5,4,3); % Matrix No.3

%% Step 2:  Maxwellian-kappa (MK) plasma matrix
% PBK type 1
Mx_mk = M_mk_mixed_pbk1(S_mk,Ns_pbk,Ns_mk,kappasz_pbk,Jpole,csnj,bx1snj,bx2snj,bx3snj,4,5,4,3); % Matrix No.4
My_mk = M_mk_mixed_pbk1(S_mk,Ns_pbk,Ns_mk,kappasz_pbk,Jpole,csnj,by1snj,by2snj,by3snj,5,5,4,3); % Matrix No.5
Mz_mk = M_mk_mixed_pbk1(S_mk,Ns_pbk,Ns_mk,kappasz_pbk,Jpole,csnj,bz1snj,bz2snj,bz3snj,6,5,4,3); % Matrix No.6

%% Step 3: Assemble matrix based on the maxied distribution of PBK and MK plasmas.
M_pbk = [Mx_pbk; My_pbk; Mz_pbk];
% disp(['size(Mx_pbk) = ',num2str(size(Mx_pbk))]);
% disp(['size(My_pbk) = ',num2str(size(My_pbk))]);
% disp(['size(Mz_pbk) = ',num2str(size(Mz_pbk))]);

M_mk = [Mx_mk; My_mk; Mz_mk];
% disp(['size(Mx_mk) = ',num2str(size(Mx_mk))]);
% disp(['size(My_mk) = ',num2str(size(My_mk))]);
% disp(['size(Mz_mk) = ',num2str(size(Mz_mk))]);

MatrixColLen_mixed = size(M_mk,2);
O = zeros(9,MatrixColLen_mixed);

M = [M_pbk; M_mk; O];

%% Step 5: The perturbed currents of first term of MK
M(end-8,end-5) = M(end-8,end-5) + bx1_mk; % the first term of dEx
M(end-8,end-4) = M(end-8,end-4) + bx2_mk; % the second term of dEy
M(end-8,end-3) = M(end-8,end-3) + bx3_mk; % the third term of dEz
%
M(end-7,end-5) = M(end-7,end-5) + by1_mk; % the first term of dEx
M(end-7,end-4) = M(end-7,end-4) + by2_mk; % the second term of dEy
M(end-7,end-3) = M(end-7,end-3) + by3_mk; % the third term of dEz
%
M(end-6,end-5) = M(end-6,end-5) + bz1_mk; % the first term of dEx
M(end-6,end-4) = M(end-6,end-4) + bz2_mk; % the second term of dEy
M(end-6,end-3) = M(end-6,end-3) + bz3_mk; % the third term of dEz

%% Step 5.1: The perturbed currents of Maxwell's equations for the PBK plasmas                 
%  Note (2026.04.02): where J_{x,y,z} is J_{x,y,z}/norm_omega
% PBK type 1
ind_Jx_pbk = getIndexOfBlkMatrix_mixed_pbk1(Ns_pbk,Ns_mk,kappasz_pbk,Jpole,1);
M(end-5,ind_Jx_pbk(end)) = M(end-5,ind_Jx_pbk(end)) - 1i/epsilon0; % for the second item of Ex

% PBK type 1
ind_Jy_pbk = getIndexOfBlkMatrix_mixed_pbk1(Ns_pbk,Ns_mk,kappasz_pbk,Jpole,2);
M(end-4,ind_Jy_pbk(end)) = M(end-4,ind_Jy_pbk(end)) - 1i/epsilon0; % for the third item of Ey

% PBK type 1
ind_Jz_pbk = getIndexOfBlkMatrix_mixed_pbk1(Ns_pbk,Ns_mk,kappasz_pbk,Jpole,3);
M(end-3,ind_Jz_pbk(end)) = M(end-3,ind_Jz_pbk(end)) - 1i/epsilon0; % for the second item of Ez

%% Step 5.2: The perturbed currents of Maxwell's equations for the Maxwellian-kappa (MK) plasmas
%  Note (2026.04.02): where J_{x,y,z} is J_{x,y,z}/norm_omega
% PBK type 1
ind_Jx_mk = getIndexOfBlkMatrix_mixed_pbk1(Ns_pbk,Ns_mk,kappasz_pbk,Jpole,4);
M(end-5,end-1) = M(end-5,end-1) + c2*kz/norm_omega; % for the first item of Ex
M(end-5,end-8) = M(end-5,end-8) - 1i/epsilon0; % for dJx1_mk
M(end-5,ind_Jx_mk(1:end-1)) = M(end-5,ind_Jx_mk(1:end-1)) - 1i/epsilon0; % for dJx2_mk

% PBK type 1
ind_Jy_mk = getIndexOfBlkMatrix_mixed_pbk1(Ns_pbk,Ns_mk,kappasz_pbk,Jpole,5);
M(end-4,end-2) = M(end-4,end-2) - c2*kz/norm_omega; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx/norm_omega; % for the second item of Ey
M(end-4,end-7) = M(end-4,end-7) - 1i/epsilon0; % for dJy1_mk
M(end-4,ind_Jy_mk(1:end-1)) = M(end-4,ind_Jy_mk(1:end-1)) - 1i/epsilon0; % for dJy2_mk

% PBK type 1
ind_Jz_mk = getIndexOfBlkMatrix_mixed_pbk1(Ns_pbk,Ns_mk,kappasz_pbk,Jpole,6);
M(end-3,end-1) = M(end-3,end-1) - c2*kx/norm_omega;  % for the first item of Ez
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0; % for dJz1_mk
M(end-3,ind_Jz_mk(1:end-1)) = M(end-3,ind_Jz_mk(1:end-1)) - 1i/epsilon0; % for dJz2_mk

%% Step 6: Maxwell's equations for the perturbed quantities of Bx, By and Bz.
M(end-2,end-4) = M(end-2,end-4) - kz/norm_omega; % for Bx
M(end-1,end-5) = M(end-1,end-5) + kz/norm_omega; % for the first item of By
M(end-1,end-3) = M(end-1,end-3) - kx/norm_omega; % for the second item of By
M(end,end-4) = M(end,end-4) + kx/norm_omega;     % for Bz

%% Solver
if(sp==0)
    % for eig()
    [V,D] = eig(double(M));
elseif(sp==1)
    % for sparse eigs();
    [V,D] = eigs(sparse(double(M)),size(M,1)); 
    % [V,D] = eigs(sparse(double(M)),20,(1+0.2i)*wps(1)); 
else
    % Convert M to a gpuArray
    M_gpu = gpuArray(double(M));
    [V,D] = eig(M_gpu);
    V = gather(V);
    D = gather(D);
end

omega = diag(D);
[~,idx]=sort(imag(omega),'descend');
% [wr,ind]=sort(real(omega),'descend');
w = omega(idx);
w(abs(w)==min(abs(w))) = NaN+1i*NaN; % remove zero solution

eigVec = V(:,idx);

end
