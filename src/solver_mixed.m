function [w,eigVec] = solver_mixed(J_opt,kx,kz,theta,B0,par,sp,EPS0)
% @Description: To solve for the wave mode roots of the oblique plasma wave 
% model using the mixed distribution in both loss-cone bi-Maxwellian and PBK 
% plasmas, the transverse wave number kx and longitudinal wave number kz 
% are specified.
% @Filename: solver_mixed.m
% @Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
% @Date: 2023-8-16
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

J = floor(J_opt);

[S,Ns,index_pbk,index_bm,kappasz,kappasx,vtsz,vtsx,Tsz,Tsx,sgms,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaParameters(B0,par);
lambdas = 0.5*kx^2*rhocs.^2;

S_pbk = sum(index_pbk);
Ns_pbk = Ns(index_pbk);
kappasz_pbk = kappasz(index_pbk);
kappasx_pbk = kappasx(index_pbk);
vtsz_pbk = vtsz(index_pbk);
vtsx_pbk = vtsx(index_pbk);
Tsz_pbk = Tsz(index_pbk);
Tsx_pbk = Tsx(index_pbk);
wps_pbk = wps(index_pbk);
wcs_pbk = wcs(index_pbk);
us0_pbk = us0(index_pbk);
rhocs_pbk = rhocs(index_pbk);
lambdaDs_pbk = lambdaDs(index_pbk);
lambdas_pbk = lambdas(index_pbk);
sgms_pbk = sgms(index_pbk); % 2025.01.22
%
S_bm = sum(index_bm);
Ns_bm = Ns(index_bm);
kappasx_bm = kappasx(index_bm);
vtsz_bm = vtsz(index_bm);
vtsx_bm = vtsx(index_bm);
Tsz_bm = Tsz(index_bm);
Tsx_bm = Tsx(index_bm);
wps_bm = wps(index_bm);
wcs_bm = wcs(index_bm);
us0_bm = us0(index_bm);
rhocs_bm = rhocs(index_bm);
lambdaDs_bm = lambdaDs(index_bm);
lambdas_bm = lambdas(index_bm);
sgms_bm = sgms(index_bm); % 2025.01.22

%distribution plasma
csn = @(s,n) pbk_csn(s,n,kz,kappasz_pbk,vtsz_pbk,wcs_pbk,us0_pbk);
b1snl = @(s,n,l) pbk_b1snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,EPS0);
b2snl = @(s,n,l) pbk_b2snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,EPS0);
b3snl = @(s,n,l) pbk_b3snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,EPS0);
b4snl = @(s,n,l) pbk_b4snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,EPS0);
b5snl = @(s,n,l) pbk_b5snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,sgms_pbk,wcs_pbk,lambdas_pbk,EPS0);
b6snl = @(s,n,l) pbk_b6snl(s,n,l,kz,kappasz_pbk,kappasx_pbk,vtsz_pbk,vtsx_pbk,sgms_pbk,lambdas_pbk,EPS0);
bx10_pbk = pbk_bx10(wps_pbk,kappasx_pbk,sgms_pbk);
bx11snl = @(s,n,l) pbk_bx11snl(s,n,l,b1snl,wps_pbk);
bx12snl = @(s,n,l) pbk_bx12snl(s,n,l,b2snl,wps_pbk);
bx21snl = @(s,n,l) pbk_bx21snl(s,n,l,b3snl,wps_pbk);
bx22snl = @(s,n,l) pbk_bx22snl(s,n,l,b4snl,wps_pbk);
bx31snl = @(s,n,l) pbk_bx31snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk);
bx32snl = @(s,n,l) pbk_bx32snl(s,n,l,theta,csn,b2snl,wps_pbk,wcs_pbk);
bx33snl = @(s,n,l) pbk_bx33snl(s,n,l,theta,b1snl,wps_pbk,wcs_pbk);
by11snl = @(s,n,l) -1*bx21snl(s,n,l);
by12snl = @(s,n,l) -1*bx22snl(s,n,l);
by20_pbk = 1i*epsilon0*sum(wps_pbk.^2); 
by21snl = @(s,n,l) pbk_by21snl(s,n,l,b5snl,wps_pbk);
by22snl = @(s,n,l) pbk_by22snl(s,n,l,b6snl,wps_pbk);
by31snl = @(s,n,l) pbk_by31snl(s,n,l,theta,csn,b3snl,b4snl,wps_pbk,wcs_pbk);
by32snl = @(s,n,l) pbk_by32snl(s,n,l,theta,csn,b4snl,wps_pbk,wcs_pbk);
by33snl = @(s,n,l) pbk_by33snl(s,n,l,theta,b3snl,wps_pbk,wcs_pbk);
bz11snl = @(s,n,l) bx31snl(s,n,l);
bz12snl = @(s,n,l) bx32snl(s,n,l);
bz13snl = @(s,n,l) bx33snl(s,n,l);
bz21snl = @(s,n,l) -1*by31snl(s,n,l);
bz22snl = @(s,n,l) -1*by32snl(s,n,l);
bz23snl = @(s,n,l) -1*by33snl(s,n,l);
bz31snl = @(s,n,l) pbk_bz31snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk);
bz32snl = @(s,n,l) pbk_bz32snl(s,n,l,theta,csn,b2snl,wps_pbk,wcs_pbk);
bz33snl = @(s,n,l) pbk_bz33snl(s,n,l,theta,csn,b1snl,b2snl,wps_pbk,wcs_pbk);
bz34snl = @(s,n,l) pbk_bz34snl(s,n,l,theta,b1snl,wps_pbk,wcs_pbk);
%
csnj = @(s,n,jj) maxwell_csnj(s,n,jj,kz,J_opt,vtsz_bm,wcs_bm,us0_bm);
b12snj = @(s,n,jj) maxwell_b12snj(s,n,jj,kz,J_opt,vtsz_bm,Tsz_bm,Tsx_bm,sgms_bm,wcs_bm,lambdas_bm,EPS0);
b34snj = @(s,n,jj) maxwell_b34snj(s,n,jj,kz,J_opt,vtsz_bm,Tsz_bm,Tsx_bm,sgms_bm,wcs_bm,lambdas_bm,EPS0);
b56snj = @(s,n,jj) maxwell_b56snj(s,n,jj,kz,J_opt,vtsz_bm,Tsz_bm,Tsx_bm,sgms_bm,wcs_bm,lambdas_bm,EPS0);
bx1_bm = maxwell_bx1(S_bm,Ns_bm,J,b12snj,csnj,wps_bm);
bx1snj = @(s,n,jj) maxwell_bx1snj(s,n,jj,b12snj,csnj,wps_bm);
bx2_bm = maxwell_bx2(S_bm,Ns_bm,J,b34snj,csnj,wps_bm);
bx2snj = @(s,n,jj) maxwell_bx2snj(s,n,jj,b34snj,csnj,wps_bm);
bx3_bm = maxwell_bx3(S_bm,Ns_bm,J,theta,b12snj,csnj,wps_bm);
bx3snj = @(s,n,jj) maxwell_bx3snj(s,n,jj,theta,b12snj,csnj,wps_bm,wcs_bm);
by1_bm = maxwell_by1(S_bm,Ns_bm,J,b34snj,csnj,wps_bm);
by1snj = @(s,n,jj) maxwell_by1snj(s,n,jj,b34snj,csnj,wps_bm);
by2_bm = maxwell_by2(S_bm,Ns_bm,J,b56snj,csnj,wps_bm);
by2snj = @(s,n,jj) maxwell_by2snj(s,n,jj,b56snj,csnj,wps_bm);
by3_bm = maxwell_by3(S_bm,Ns_bm,J,theta,b34snj,csnj,wps_bm);
by3snj = @(s,n,jj) maxwell_by3snj(s,n,jj,theta,b34snj,csnj,wps_bm,wcs_bm);
bz1_bm = maxwell_bz1(S_bm,Ns_bm,J,theta,b12snj,csnj,wps_bm);
bz1snj = @(s,n,jj) maxwell_bz1snj(s,n,jj,theta,b12snj,csnj,wps_bm,wcs_bm);
bz2_bm = maxwell_bz2(S_bm,Ns_bm,J,theta,b34snj,csnj,wps_bm);
bz2snj = @(s,n,jj) maxwell_bz2snj(s,n,jj,theta,b34snj,csnj,wps_bm,wcs_bm);
bz3_bm = maxwell_bz3(S_bm,Ns_bm,J,theta,b12snj,csnj,wps_bm);
bz3snj = @(s,n,jj) maxwell_bz3snj(s,n,jj,theta,b12snj,csnj,wps_bm,wcs_bm);

% matrix
Mx_pbk = Mxy_pbk_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csn, ...
                       bx11snl,bx21snl,bx31snl,bx12snl,bx22snl,bx32snl, ...
                       bx33snl,bx10_pbk,1,5,4,3,5); % Matrix No.1
My_pbk = Mxy_pbk_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csn, ...
                       by11snl,by21snl,by31snl,by12snl,by22snl,by32snl, ...
                       by33snl,by20_pbk,2,5,4,3,4); % Matrix No.2
Mz_pbk = Mz_pbk_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csn, ...
                      bz11snl,bz21snl,bz31snl,bz13snl,bz23snl,bz33snl, ...
                      bz12snl,bz22snl,bz32snl,bz34snl,by20_pbk,3,5,4,3); % Matrix No.3

Mx_bm = M_mixed_maxwell(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csnj,bx1snj,bx2snj,bx3snj,4,5,4,3); % Matrix No.4
My_bm = M_mixed_maxwell(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csnj,by1snj,by2snj,by3snj,5,5,4,3); % Matrix No.5
Mz_bm = M_mixed_maxwell(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csnj,bz1snj,bz2snj,bz3snj,6,5,4,3); % Matrix No.6
M_pbk = [Mx_pbk; My_pbk; Mz_pbk];
M_bm = [Mx_bm; My_bm; Mz_bm];

MatrixColLen_mixed = size(M_bm,2);
O = zeros(9,MatrixColLen_mixed);
M = [M_pbk; M_bm; O];
M(end-8,end-5) = M(end-8,end-5) + bx1_bm;
M(end-8,end-4) = M(end-8,end-4) + bx2_bm; 
M(end-8,end-3) = M(end-8,end-3) + bx3_bm;
M(end-7,end-5) = M(end-7,end-5) + by1_bm;
M(end-7,end-4) = M(end-7,end-4) + by2_bm;
M(end-7,end-3) = M(end-7,end-3) + by3_bm;
M(end-6,end-5) = M(end-6,end-5) + bz1_bm;
M(end-6,end-4) = M(end-6,end-4) + bz2_bm;
M(end-6,end-3) = M(end-6,end-3) + bz3_bm;                  
idx_Jx_pbk = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,1);
M(end-5,idx_Jx_pbk(end)) = M(end-5,idx_Jx_pbk(end)) - 1i/epsilon0; % for the second item of Ex
idx_Jy_pbk = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,2);
M(end-4,idx_Jy_pbk(end)) = M(end-4,idx_Jy_pbk(end)) - 1i/epsilon0; % for the third item of Ey
idx_Jz_pbk = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,3);
M(end-3,idx_Jz_pbk(end)) = M(end-3,idx_Jz_pbk(end)) - 1i/epsilon0; % for the second item of Ez
idx_Jx_bm = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,4);
M(end-5,end-1) = M(end-5,end-1) + c2*kz; % for the first item of Ex
M(end-5,end-8) = M(end-5,end-8) - 1i/epsilon0; % for dJx1_bm
M(end-5,idx_Jx_bm(1:end-1)) = M(end-5,idx_Jx_bm(1:end-1)) - 1i/epsilon0; % for dJx2_bm
idx_Jy_bm = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,5);
M(end-4,end-2) = M(end-4,end-2) - c2*kz; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx; % for the second item of Ey
M(end-4,end-7) = M(end-4,end-7) - 1i/epsilon0; % for dJy1_bm
M(end-4,idx_Jy_bm(1:end-1)) = M(end-4,idx_Jy_bm(1:end-1)) - 1i/epsilon0; % for dJy2_bm
idx_Jz_bm = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,6);
M(end-3,end-1) = M(end-3,end-1) - c2*kx;  % for the first item of Ez
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0; % for dJz1_bm
M(end-3,idx_Jz_bm(1:end-1)) = M(end-3,idx_Jz_bm(1:end-1)) - 1i/epsilon0; % for dJz2_bm
M(end-2,end-4) = M(end-2,end-4) - kz;
M(end-1,end-5) = M(end-1,end-5) + kz; 
M(end-1,end-3) = M(end-1,end-3) - kx;
M(end,end-4) = M(end,end-4) + kx;

%% Solver
if(sp==0)
    % for eig()
    [V,D] = eig(double(M));
elseif(sp==1)
    % for sparse eigs();
    [V,D] = eigs(sparse(double(M)),size(M,1)); 
else
    % Convert M to a gpuArray
    M_gpu = gpuArray(double(M));
    [V,D] = eig(M_gpu);
    V = gather(V);
    D = gather(D);
end

omega = diag(D);
[~,idx]=sort(imag(omega),'descend');
w = omega(idx);
w(abs(w)==min(abs(w))) = NaN+1i*NaN; 

eigVec = V(:,idx);

end
