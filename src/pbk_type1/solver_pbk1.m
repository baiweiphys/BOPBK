function [w,eigVec] = solver_pbk1(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0)
% SOLVER_PBK To solve for the wave mode roots in the PBK solver. 
% @Description: To solve for the wave mode roots of the oblique plasma wave 
% model using the loss-cone PBK (type 1) distribution, the transverse wave number kx
% and longitudinal wave number kz are specified.
% @Filename: solver_pbk1.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-9-24
% @Modification: 
%   Bai Wei, 2026.04.01: add norm_omega.
%   Bai Wei, 2026.04.05: Update solver to PBK type1.

params_with_unit;

% Input parameters
S = plasmaParams.S;
Ns = plasmaParams.Ns;
% J_opt = plasmaParams.J_opt;
% index_pbk = plasmaParams.index_pbk;
% index_mk = plasmaParams.index_mk;
kappasz = plasmaParams.kappasz;
kappasx = plasmaParams.kappasx;
vtsz = plasmaParams.vtsz;
vtsx = plasmaParams.vtsx;
% Tsz = plasmaParams.Tsz;
% Tsx = plasmaParams.Tsx;
sgms = plasmaParams.sgms;
wps = plasmaParams.wps;
wcs = plasmaParams.wcs;
us0 = plasmaParams.us0;
rhocs = plasmaParams.rhocs;
% is_pbk_vts = plasmaParams.is_pbk_vts;
% lambdaDs = plasmaParams.lambdaDs;

lambdas = 0.5*kx.^2*rhocs.^2; % for argument of the modified bessel function
% lambdas(abs(lambdas)<eps) = eps;  % 2024.0916, to avoid singular when k_perp=0

% common/PBK type 1 or 2
csn = @(s,n) pbk_csn(s,n,kz,kappasz,vtsz,wcs,us0,norm_omega);

% PBK type 1
b1snl = @(s,n,l) pbk1_b1snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,norm_omega,EPS0);
b2snl = @(s,n,l) pbk1_b2snl(s,n,l,kz,kappasz,kappasx,vtsz,vtsx,sgms,lambdas,norm_omega,EPS0);
b3snl = @(s,n,l) pbk1_b3snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,norm_omega,EPS0);
b4snl = @(s,n,l) pbk1_b4snl(s,n,l,kz,kappasz,kappasx,vtsz,vtsx,sgms,lambdas,norm_omega,EPS0);
b5snl = @(s,n,l) pbk1_b5snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,norm_omega,EPS0);
b6snl = @(s,n,l) pbk1_b6snl(s,n,l,kz,kappasz,kappasx,vtsz,vtsx,sgms,lambdas,norm_omega,EPS0);

% for x-component (common/PBK type 1 or 2)
bx10_pbk = pbk_bx10(wps,kappasx,sgms,norm_omega);
% bx30_pbk = pbk_bx30(theta,wps);
bx11snl = @(s,n,l) pbk_bx11snl(s,n,l,b1snl,wps,norm_omega);
bx12snl = @(s,n,l) pbk_bx12snl(s,n,l,b2snl,wps,norm_omega);
bx21snl = @(s,n,l) pbk_bx21snl(s,n,l,b3snl,wps,norm_omega);
bx22snl = @(s,n,l) pbk_bx22snl(s,n,l,b4snl,wps,norm_omega);
bx31snl = @(s,n,l) pbk_bx31snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs,norm_omega);
bx32snl = @(s,n,l) pbk_bx32snl(s,n,l,theta,csn,b2snl,wps,wcs,norm_omega);
bx33snl = @(s,n,l) pbk_bx33snl(s,n,l,theta,b1snl,wps,wcs,norm_omega);

% for y-component (common/PBK type 1 or 2)
by11snl = @(s,n,l) -1*bx21snl(s,n,l);
by12snl = @(s,n,l) -1*bx22snl(s,n,l);
by20_pbk = 1i*epsilon0*sum(wps.^2)/norm_omega^2;
by21snl = @(s,n,l) pbk_by21snl(s,n,l,b5snl,wps,norm_omega);
by22snl = @(s,n,l) pbk_by22snl(s,n,l,b6snl,wps,norm_omega);
by31snl = @(s,n,l) pbk_by31snl(s,n,l,theta,csn,b3snl,b4snl,wps,wcs,norm_omega);
by32snl = @(s,n,l) pbk_by32snl(s,n,l,theta,csn,b4snl,wps,wcs,norm_omega);
by33snl = @(s,n,l) pbk_by33snl(s,n,l,theta,b3snl,wps,wcs,norm_omega);

% for z-component (common/PBK type 1 or 2)
% bz10_pbk = bx30_pbk;
bz11snl = @(s,n,l) bx31snl(s,n,l);
bz12snl = @(s,n,l) bx32snl(s,n,l);
bz13snl = @(s,n,l) bx33snl(s,n,l);
bz21snl = @(s,n,l) -1*by31snl(s,n,l);
bz22snl = @(s,n,l) -1*by32snl(s,n,l);
bz23snl = @(s,n,l) -1*by33snl(s,n,l);
% bz30_pbk = pbk_bz30(theta,wps);
bz31snl = @(s,n,l) pbk_bz31snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs,norm_omega);
bz32snl = @(s,n,l) pbk_bz32snl(s,n,l,theta,csn,b2snl,wps,wcs,norm_omega);
bz33snl = @(s,n,l) pbk_bz33snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs,norm_omega);
bz34snl = @(s,n,l) pbk_bz34snl(s,n,l,theta,b1snl,wps,wcs,norm_omega);


% Step 1
% PBK type 1
Mx = Mxy_pbk1(S,Ns,kappasz,csn,bx11snl,bx21snl,bx31snl, ...
            bx12snl,bx22snl,bx32snl,bx33snl,bx10_pbk,1,5,4,3,5); % Matrix No.1

% PBK type 1
My = Mxy_pbk1(S,Ns,kappasz,csn,by11snl,by21snl,by31snl, ...
             by12snl,by22snl,by32snl,by33snl,by20_pbk,2,5,4,3,4); % Matrix No.2
 
% PBK type 1
Mz = Mz_pbk1(S,Ns,kappasz,csn,bz11snl,bz21snl,bz31snl, ...
            bz13snl,bz23snl,bz33snl,bz12snl,bz22snl,bz32snl,bz34snl,by20_pbk,3,5,4,3); % Matrix No.3


% Step 2
MatrixColLen = size(Mz,2);
O = zeros(6,MatrixColLen);

M = [Mx; My; Mz; O];

%% Step 3: Maxwell's equations for Ex
%  Note (2026.04.02): where J_{x,y,z} is J_{x,y,z}/norm_omega
% PBK type 1
ind_Jx_pbk = getIndexOfBlkMatrix_pbk1(Ns,kappasz,1);
M(end-5,end-1) = M(end-5,end-1) + c2*kz/norm_omega; % for the first item of Ex
M(end-5,ind_Jx_pbk(end)) = M(end-5,ind_Jx_pbk(end)) - 1i/epsilon0; % for the second item of Ex

% PBK type 1
ind_Jy_pbk = getIndexOfBlkMatrix_pbk1(Ns,kappasz,2);
M(end-4,end-2) = M(end-4,end-2) - c2*kz/norm_omega; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx/norm_omega; % for the second item of Ey
M(end-4,ind_Jy_pbk(end)) = M(end-4,ind_Jy_pbk(end)) - 1i/epsilon0; % for the third item of Ey

% PBK type 1
ind_Jz_pbk = getIndexOfBlkMatrix_pbk1(Ns,kappasz,3);
M(end-3,end-1) = M(end-3,end-1) - c2*kx/norm_omega;  % for the first item of Ez
M(end-3,ind_Jz_pbk(end)) = M(end-3,ind_Jz_pbk(end)) - 1i/epsilon0; % for the second item of Ez

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
    % [V,D] = eigs(sparse(double(M)),100,'largestimag'); 
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
