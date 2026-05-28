function  [w,eigVec] = solver_mk1(kx,kz,theta,plasmaParams,sp,norm_omega,EPS0)
% SOLVER_MK1 To solve for the wave mode roots in the MK solver. 
% @Description: To solve for the wave mode roots of the oblique plasma wave 
% model using the loss-cone Maxwellian-kappa (MK, PBK type 1) distribution, the transverse wave 
% number kx and longitudinal wave number kz are specified.
% @Filename: solver_mk1.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-9-03
% @Modification: 
%   Bai Wei, 2026.04.02: add norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM structure).
%                        Uses BM properties only if (kappaxs > kappaxs_th); 
%                        otherwise uses MK.
%   Bai Wei, 2026.04.05: Update solver to MK (PBK type1)

params_with_unit;

% Input parameters
S = plasmaParams.S;
Ns = plasmaParams.Ns;
J_opt = plasmaParams.J_opt;
% index_pbk = plasmaParams.index_pbk;
% index_mk = plasmaParams.index_mk;
% kappasz = plasmaParams.kappasz;
kappasx = plasmaParams.kappasx;
kappasx_th = plasmaParams.kappasx_threshold;
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

Jpole = floor(J_opt);  % Use the floor function to get the integer part

%%
lambdas = 0.5*kx.^2*rhocs.^2; % for argument of the modified bessel function
lambdas(abs(lambdas)<1e-50) = 1e-50;  % 2024.0916, to avoid singular when k_perp=0

% common/PBK type 1 or 2
csnj = @(s,n,jj) mk_csnj(s,n,jj,kz,J_opt,vtsz,wcs,us0,norm_omega);

% PBK type 1
b12snj = @(s,n,jj) mk1_b12snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,norm_omega,EPS0);
b34snj = @(s,n,jj) mk1_b34snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,norm_omega,EPS0);
b56snj = @(s,n,jj) mk1_b56snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,norm_omega,EPS0);

% common/PBK type 1 or 2
bx1_mk = mk_bx1(S,Ns,Jpole,b12snj,csnj,wps,norm_omega);
bx1snj = @(s,n,jj) mk_bx1snj(s,n,jj,b12snj,csnj,wps,norm_omega);
bx2_mk = mk_bx2(S,Ns,Jpole,b34snj,csnj,wps,norm_omega);
bx2snj = @(s,n,jj) mk_bx2snj(s,n,jj,b34snj,csnj,wps,norm_omega);
bx3_mk = mk_bx3(S,Ns,Jpole,theta,b12snj,csnj,wps,norm_omega);
bx3snj = @(s,n,jj) mk_bx3snj(s,n,jj,theta,b12snj,csnj,wps,wcs,norm_omega);

% common/PBK type 1 or 2
by1_mk = mk_by1(S,Ns,Jpole,b34snj,csnj,wps,norm_omega);
by1snj = @(s,n,jj) mk_by1snj(s,n,jj,b34snj,csnj,wps,norm_omega);
by2_mk = mk_by2(S,Ns,Jpole,b56snj,csnj,wps,norm_omega);
by2snj = @(s,n,jj) mk_by2snj(s,n,jj,b56snj,csnj,wps,norm_omega);
by3_mk = mk_by3(S,Ns,Jpole,theta,b34snj,csnj,wps,norm_omega);
by3snj = @(s,n,jj) mk_by3snj(s,n,jj,theta,b34snj,csnj,wps,wcs,norm_omega);

% common/PBK type 1 or 2
bz1_mk = mk_bz1(S,Ns,Jpole,theta,b12snj,csnj,wps,norm_omega);
bz1snj = @(s,n,jj) mk_bz1snj(s,n,jj,theta,b12snj,csnj,wps,wcs,norm_omega);
bz2_mk = mk_bz2(S,Ns,Jpole,theta,b34snj,csnj,wps,norm_omega);
bz2snj = @(s,n,jj) mk_bz2snj(s,n,jj,theta,b34snj,csnj,wps,wcs,norm_omega);
bz3_mk = mk_bz3(S,Ns,Jpole,theta,b12snj,csnj,wps,norm_omega);
bz3snj = @(s,n,jj) mk_bz3snj(s,n,jj,theta,b12snj,csnj,wps,wcs,norm_omega);

% Step 1
% for x-component (common/PBK type 1 or 2)
Mx_mk = M_mk(S,Ns,Jpole,csnj,bx1snj,bx2snj,bx3snj,1,5,4,3); % Matrix No.1

% for y-component (common/PBK type 1 or 2)
My_mk = M_mk(S,Ns,Jpole,csnj,by1snj,by2snj,by3snj,2,5,4,3); % Matrix No.2

% for z-component (common/PBK type 1 or 2)
Mz_mk = M_mk(S,Ns,Jpole,csnj,bz1snj,bz2snj,bz3snj,3,5,4,3); % Matrix No.3

% Step 2
MatrixColLen = size(Mx_mk,2);
O = zeros(9,MatrixColLen);

M = [Mx_mk; My_mk; Mz_mk; O];

%% Step 3: The perturbed currents of first term of Maxwellian-kappa (MK)
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

%% Step 4: Maxwell's equation
%  Note (2026.04.02): where J_{x,y,z} is J_{x,y,z}/norm_omega
% common/PBK type 1 or 2
ind_Jx_mk = getIndexOfBlkMatrix_mk(Ns,Jpole,1);
M(end-5,end-1) = M(end-5,end-1) + c2*kz/norm_omega; % for the first item of Ex
M(end-5,end-8) = M(end-5,end-8) - 1i/epsilon0; % for dJx1
M(end-5,ind_Jx_mk(1:end-1)) = M(end-5,ind_Jx_mk(1:end-1)) - 1i/epsilon0; % for dJx2
%
% common/PBK type 1 or 2
ind_Jy_mk = getIndexOfBlkMatrix_mk(Ns,Jpole,2); 
M(end-4,end-2) = M(end-4,end-2) - c2*kz/norm_omega; % for the first item of Ey
M(end-4,end) = M(end-4,end) + c2*kx/norm_omega; % for the second item of Ey
M(end-4,end-7) = M(end-4,end-7) - 1i/epsilon0; % for dJy1
M(end-4,ind_Jy_mk(1:end-1)) = M(end-4,ind_Jy_mk(1:end-1)) - 1i/epsilon0; % for dJy2
%
% common/PBK type 1 or 2
ind_Jz_mk = getIndexOfBlkMatrix_mk(Ns,Jpole,3);
M(end-3,end-1) = M(end-3,end-1) - c2*kx/norm_omega;  % for the first item of Ez
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0; % for dJz1
M(end-3,ind_Jz_mk(1:end-1)) = M(end-3,ind_Jz_mk(1:end-1)) - 1i/epsilon0; % for dJz2
%
M(end-2,end-4) = M(end-2,end-4) - kz/norm_omega; % for Bx
M(end-1,end-5) = M(end-1,end-5) + kz/norm_omega; % for the first item of By
M(end-1,end-3) = M(end-1,end-3) - kx/norm_omega; % for the second item of By
M(end,end-4) = M(end,end-4) + kx/norm_omega;     % for Bz


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
