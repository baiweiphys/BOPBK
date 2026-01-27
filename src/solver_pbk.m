function [w,eigVec] = solver_pbk(kx,kz,theta,plasmaParams,sp,EPS0)
% @Description: To solve for the wave mode roots of the oblique plasma wave 
% model using the loss-cone PBK distribution, the transverse wave number kx
% and longitudinal wave number kz are specified.
% @Filename: solver_pbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-9-24
% @LastEditors: Bai Wei
% @LastEditTime: 2026-01-24


params_with_unit;

% Input parameters
S = plasmaParams.S;
Ns = plasmaParams.Ns;
% J_opt = plasmaParams.J_opt;
% index_pbk = plasmaParams.index_pbk;
% index_bm = plasmaParams.index_bm;
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
% lambdaDs = plasmaParams.lambdaDs;

lambdas = 0.5*kx.^2*rhocs.^2; % for argument of the modified bessel function
% lambdas(abs(lambdas)<eps) = eps;  % 2024.0916, to avoid singular when k_perp=0
csn = @(s,n) pbk_csn(s,n,kz,kappasz,vtsz,wcs,us0);
b1snl = @(s,n,l) pbk_b1snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,EPS0);
b2snl = @(s,n,l) pbk_b2snl(s,n,l,kz,kappasz,kappasx,vtsz,vtsx,sgms,lambdas,EPS0);
b3snl = @(s,n,l) pbk_b3snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,EPS0);
b4snl = @(s,n,l) pbk_b4snl(s,n,l,kz,kappasz,kappasx,vtsz,vtsx,sgms,lambdas,EPS0);
b5snl = @(s,n,l) pbk_b5snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,EPS0);
b6snl = @(s,n,l) pbk_b6snl(s,n,l,kz,kappasz,kappasx,vtsz,vtsx,sgms,lambdas,EPS0);

% for x-component
bx10_pbk = pbk_bx10(wps,kappasx,sgms);
% bx30_pbk = pbk_bx30(theta,wps);
bx11snl = @(s,n,l) pbk_bx11snl(s,n,l,b1snl,wps);
bx12snl = @(s,n,l) pbk_bx12snl(s,n,l,b2snl,wps);
bx21snl = @(s,n,l) pbk_bx21snl(s,n,l,b3snl,wps);
bx22snl = @(s,n,l) pbk_bx22snl(s,n,l,b4snl,wps);
bx31snl = @(s,n,l) pbk_bx31snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs);
bx32snl = @(s,n,l) pbk_bx32snl(s,n,l,theta,csn,b2snl,wps,wcs);
bx33snl = @(s,n,l) pbk_bx33snl(s,n,l,theta,b1snl,wps,wcs);

% for y-component
by11snl = @(s,n,l) -1*bx21snl(s,n,l);
by12snl = @(s,n,l) -1*bx22snl(s,n,l);
by20_pbk = 1i*epsilon0*sum(wps.^2);
by21snl = @(s,n,l) pbk_by21snl(s,n,l,b5snl,wps);
by22snl = @(s,n,l) pbk_by22snl(s,n,l,b6snl,wps);
by31snl = @(s,n,l) pbk_by31snl(s,n,l,theta,csn,b3snl,b4snl,wps,wcs);
by32snl = @(s,n,l) pbk_by32snl(s,n,l,theta,csn,b4snl,wps,wcs);
by33snl = @(s,n,l) pbk_by33snl(s,n,l,theta,b3snl,wps,wcs);

% for z-component
% bz10_pbk = bx30_pbk;
bz11snl = @(s,n,l) bx31snl(s,n,l);
bz12snl = @(s,n,l) bx32snl(s,n,l);
bz13snl = @(s,n,l) bx33snl(s,n,l);
bz21snl = @(s,n,l) -1*by31snl(s,n,l);
bz22snl = @(s,n,l) -1*by32snl(s,n,l);
bz23snl = @(s,n,l) -1*by33snl(s,n,l);
% bz30_pbk = pbk_bz30(theta,wps);
bz31snl = @(s,n,l) pbk_bz31snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs);
bz32snl = @(s,n,l) pbk_bz32snl(s,n,l,theta,csn,b2snl,wps,wcs);
bz33snl = @(s,n,l) pbk_bz33snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs);
bz34snl = @(s,n,l) pbk_bz34snl(s,n,l,theta,b1snl,wps,wcs);


% Step 1
Mx = Mxy_pbk(S,Ns,kappasz,csn,bx11snl,bx21snl,bx31snl, ...
             bx12snl,bx22snl,bx32snl,bx33snl,bx10_pbk,1,5,4,3,5); % Matrix No.1

My = Mxy_pbk(S,Ns,kappasz,csn,by11snl,by21snl,by31snl, ...
             by12snl,by22snl,by32snl,by33snl,by20_pbk,2,5,4,3,4); % Matrix No.2
 
Mz = Mz_pbk(S,Ns,kappasz,csn,bz11snl,bz21snl,bz31snl, ...
            bz13snl,bz23snl,bz33snl,bz12snl,bz22snl,bz32snl,bz34snl,by20_pbk,3,5,4,3); % Matrix No.3


% Step 2
MatrixColLen = size(Mz,2);
O = zeros(6,MatrixColLen);

M = [Mx; My; Mz; O];

%% Step 3: Maxwell's equations for Ex
idx_Jx_pbk = getIndexOfBlkMatrix_pbk(Ns,kappasz,1);
M(end-5,end-1) = M(end-5,end-1) + c2*kz; % for the first item of Ex
M(end-5,idx_Jx_pbk(end)) = M(end-5,idx_Jx_pbk(end)) - 1i/epsilon0; % for the second item of Ex

idx_Jy_pbk = getIndexOfBlkMatrix_pbk(Ns,kappasz,2);
M(end-4,end-2) = M(end-4,end-2) - c2*kz; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx; % for the second item of Ey
M(end-4,idx_Jy_pbk(end)) = M(end-4,idx_Jy_pbk(end)) - 1i/epsilon0; % for the third item of Ey

idx_Jz_pbk = getIndexOfBlkMatrix_pbk(Ns,kappasz,3);
M(end-3,end-1) = M(end-3,end-1) - c2*kx;  % for the first item of Ez
M(end-3,idx_Jz_pbk(end)) = M(end-3,idx_Jz_pbk(end)) - 1i/epsilon0; % for the second item of Ez

M(end-2,end-4) = M(end-2,end-4) - kz; % for Bx
M(end-1,end-5) = M(end-1,end-5) + kz; % for the first item of By
M(end-1,end-3) = M(end-1,end-3) - kx; % for the second item of By
M(end,end-4) = M(end,end-4) + kx;     % for Bz

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