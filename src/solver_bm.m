function  [w,eigVec] = solver_bm(kx,kz,theta,plasmaParams,sp,EPS0)
% @Description: To solve for the wave mode roots of the oblique plasma wave 
% model using the loss-cone bi-Maxwellian distribution, the transverse wave 
% number kx and longitudinal wave number kz are specified.
% @Filename: solver_bm.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-9-03
% @LastEditors: Bai Wei
% @LastEditTime: 2026-01-24

params_with_unit;

% Input parameters
S = plasmaParams.S;
Ns = plasmaParams.Ns;
J_opt = plasmaParams.J_opt;
% index_pbk = plasmaParams.index_pbk;
% index_bm = plasmaParams.index_bm;
% kappasz = plasmaParams.kappasz;
% kappasx = plasmaParams.kappasx;
vtsz = plasmaParams.vtsz;
% vtsx = plasmaParams.vtsx;
Tsz = plasmaParams.Tsz;
Tsx = plasmaParams.Tsx;
sgms = plasmaParams.sgms;
wps = plasmaParams.wps;
wcs = plasmaParams.wcs;
us0 = plasmaParams.us0;
rhocs = plasmaParams.rhocs;
% lambdaDs = plasmaParams.lambdaDs;

J = floor(J_opt);  % Use the floor function to get the integer part

%%
lambdas = 0.5*kx.^2*rhocs.^2; % for argument of the modified bessel function
lambdas(abs(lambdas)<1e-50) = 1e-50;  % 2024.0916, to avoid singular when k_perp=0

csnj = @(s,n,jj) maxwell_csnj(s,n,jj,kz,J_opt,vtsz,wcs,us0);

b12snj = @(s,n,jj) maxwell_b12snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0);
b34snj = @(s,n,jj) maxwell_b34snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0);
b56snj = @(s,n,jj) maxwell_b56snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0);

%
bx1_bm = maxwell_bx1(S,Ns,J,b12snj,csnj,wps);
bx1snj = @(s,n,jj) maxwell_bx1snj(s,n,jj,b12snj,csnj,wps);
bx2_bm = maxwell_bx2(S,Ns,J,b34snj,csnj,wps);
bx2snj = @(s,n,jj) maxwell_bx2snj(s,n,jj,b34snj,csnj,wps);
bx3_bm = maxwell_bx3(S,Ns,J,theta,b12snj,csnj,wps);
bx3snj = @(s,n,jj) maxwell_bx3snj(s,n,jj,theta,b12snj,csnj,wps,wcs);
%
by1_bm = maxwell_by1(S,Ns,J,b34snj,csnj,wps);
by1snj = @(s,n,jj) maxwell_by1snj(s,n,jj,b34snj,csnj,wps);
by2_bm = maxwell_by2(S,Ns,J,b56snj,csnj,wps);
by2snj = @(s,n,jj) maxwell_by2snj(s,n,jj,b56snj,csnj,wps);
by3_bm = maxwell_by3(S,Ns,J,theta,b34snj,csnj,wps);
by3snj = @(s,n,jj) maxwell_by3snj(s,n,jj,theta,b34snj,csnj,wps,wcs);
%
bz1_bm = maxwell_bz1(S,Ns,J,theta,b12snj,csnj,wps);
bz1snj = @(s,n,jj) maxwell_bz1snj(s,n,jj,theta,b12snj,csnj,wps,wcs);
bz2_bm = maxwell_bz2(S,Ns,J,theta,b34snj,csnj,wps);
bz2snj = @(s,n,jj) maxwell_bz2snj(s,n,jj,theta,b34snj,csnj,wps,wcs);
bz3_bm = maxwell_bz3(S,Ns,J,theta,b12snj,csnj,wps);
bz3snj = @(s,n,jj) maxwell_bz3snj(s,n,jj,theta,b12snj,csnj,wps,wcs);

% Step 1
% for x-component
Mx_bm = M_maxwell(S,Ns,J,csnj,bx1snj,bx2snj,bx3snj,1,5,4,3); % Matrix No.1

% for y-component
My_bm = M_maxwell(S,Ns,J,csnj,by1snj,by2snj,by3snj,2,5,4,3); % Matrix No.2

% for z-component
Mz_bm = M_maxwell(S,Ns,J,csnj,bz1snj,bz2snj,bz3snj,3,5,4,3); % Matrix No.3

% Step 2
MatrixColLen = size(Mx_bm,2);
O = zeros(9,MatrixColLen);

M = [Mx_bm; My_bm; Mz_bm; O];

%% Step 3: The perturbed currents of first term of BM

%
M(end-8,end-5) = M(end-8,end-5) + bx1_bm; % the first term of dEx
M(end-8,end-4) = M(end-8,end-4) + bx2_bm; % the second term of dEy
M(end-8,end-3) = M(end-8,end-3) + bx3_bm; % the third term of dEz
%
M(end-7,end-5) = M(end-7,end-5) + by1_bm; % the first term of dEx
M(end-7,end-4) = M(end-7,end-4) + by2_bm; % the second term of dEy
M(end-7,end-3) = M(end-7,end-3) + by3_bm; % the third term of dEz

%
M(end-6,end-5) = M(end-6,end-5) + bz1_bm; % the first term of dEx
M(end-6,end-4) = M(end-6,end-4) + bz2_bm; % the second term of dEy
M(end-6,end-3) = M(end-6,end-3) + bz3_bm; % the third term of dEz

%% Step 4: Maxwell's equation
%
idx_Jx_bm = getIndexOfBlkMatrix_maxwell(Ns,J,1);
M(end-5,end-1) = M(end-5,end-1) + c2*kz; % for the first item of Ex
M(end-5,end-8) = M(end-5,end-8) - 1i/epsilon0; % for dJx1
M(end-5,idx_Jx_bm(1:end-1)) = M(end-5,idx_Jx_bm(1:end-1)) - 1i/epsilon0; % for dJx2
%
idx_Jy_bm = getIndexOfBlkMatrix_maxwell(Ns,J,2);
M(end-4,end-2) = M(end-4,end-2) - c2*kz; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx; % for the second item of Ey
M(end-4,end-7) = M(end-4,end-7) - 1i/epsilon0; % for dJy1
M(end-4,idx_Jy_bm(1:end-1)) = M(end-4,idx_Jy_bm(1:end-1)) - 1i/epsilon0; % for dJy2
%
idx_Jz_bm = getIndexOfBlkMatrix_maxwell(Ns,J,3);
M(end-3,end-1) = M(end-3,end-1) - c2*kx;  % for the first item of Ez
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0; % for dJz1
M(end-3,idx_Jz_bm(1:end-1)) = M(end-3,idx_Jz_bm(1:end-1)) - 1i/epsilon0; % for dJz2
%
M(end-2,end-4) = M(end-2,end-4) - kz; % for Bx
M(end-1,end-5) = M(end-1,end-5) + kz; % for the first item of By
M(end-1,end-3) = M(end-1,end-3) - kx; % for the second item of By
M(end,end-4) = M(end,end-4) + kx;     % for Bz


%% %%%%%%%% Solver
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
% [wr,ind]=sort(real(omega),'descend');
w = omega(idx);
w(abs(w)==min(abs(w))) = NaN+1i*NaN; % remove zero solution

eigVec = V(:,idx);

end