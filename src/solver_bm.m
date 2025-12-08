function  [w,eigVec] = solver_bm(J_opt,kx,kz,theta,B0,par,sp,EPS0)
% @Description: To solve for the wave mode roots of the oblique plasma wave 
% model using the loss-cone bi-Maxwellian distribution, the transverse wave 
% number kx and longitudinal wave number kz are specified.
% @Filename: solver_bm.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-9-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-06

params_with_unit;

J = floor(J_opt);

[S,Ns,index_pbk,index_bm,kappasz,kappasx,vtsz,vtsx,Tsz,Tsx,sgms,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaParameters(B0,par); 

lambdas = 0.5*kx.^2*rhocs.^2;
lambdas(abs(lambdas)<1e-50) = 1e-50;  % 2024.0916, to avoid singular when k_perp=0

csnj = @(s,n,jj) maxwell_csnj(s,n,jj,kz,J_opt,vtsz,wcs,us0);
b12snj = @(s,n,jj) maxwell_b12snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0);
b34snj = @(s,n,jj) maxwell_b34snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0);
b56snj = @(s,n,jj) maxwell_b56snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0);

bx1_bm = maxwell_bx1(S,Ns,J,b12snj,csnj,wps);
bx1snj = @(s,n,jj) maxwell_bx1snj(s,n,jj,b12snj,csnj,wps);
bx2_bm = maxwell_bx2(S,Ns,J,b34snj,csnj,wps);
bx2snj = @(s,n,jj) maxwell_bx2snj(s,n,jj,b34snj,csnj,wps);
bx3_bm = maxwell_bx3(S,Ns,J,theta,b12snj,csnj,wps);
bx3snj = @(s,n,jj) maxwell_bx3snj(s,n,jj,theta,b12snj,csnj,wps,wcs);
by1_bm = maxwell_by1(S,Ns,J,b34snj,csnj,wps);
by1snj = @(s,n,jj) maxwell_by1snj(s,n,jj,b34snj,csnj,wps);
by2_bm = maxwell_by2(S,Ns,J,b56snj,csnj,wps);
by2snj = @(s,n,jj) maxwell_by2snj(s,n,jj,b56snj,csnj,wps);
by3_bm = maxwell_by3(S,Ns,J,theta,b34snj,csnj,wps);
by3snj = @(s,n,jj) maxwell_by3snj(s,n,jj,theta,b34snj,csnj,wps,wcs);
bz1_bm = maxwell_bz1(S,Ns,J,theta,b12snj,csnj,wps);
bz1snj = @(s,n,jj) maxwell_bz1snj(s,n,jj,theta,b12snj,csnj,wps,wcs);
bz2_bm = maxwell_bz2(S,Ns,J,theta,b34snj,csnj,wps);
bz2snj = @(s,n,jj) maxwell_bz2snj(s,n,jj,theta,b34snj,csnj,wps,wcs);
bz3_bm = maxwell_bz3(S,Ns,J,theta,b12snj,csnj,wps);
bz3snj = @(s,n,jj) maxwell_bz3snj(s,n,jj,theta,b12snj,csnj,wps,wcs);

% matrix
by20 = 1i*epsilon0*sum(wps.^2);
Mx_bm = M_maxwell(S,Ns,J,csnj,bx1snj,bx2snj,bx3snj,1,5,4,3); % Matrix No.1
My_bm = M_maxwell(S,Ns,J,csnj,by1snj,by2snj,by3snj,2,5,4,3); % Matrix No.2
Mz_bm = M_maxwell(S,Ns,J,csnj,bz1snj,bz2snj,bz3snj,3,5,4,3); % Matrix No.3

MatrixColLen = size(Mx_bm,2);
O = zeros(9,MatrixColLen);

M = [Mx_bm; My_bm; Mz_bm; O];

M(end-8,end-5) = M(end-8,end-5) + bx1_bm; 
M(end-8,end-4) = M(end-8,end-4) + bx2_bm;
M(end-8,end-3) = M(end-8,end-3) + bx3_bm; 
M(end-7,end-5) = M(end-7,end-5) + by1_bm; 
M(end-7,end-4) = M(end-7,end-4) + by2_bm;
M(end-7,end-3) = M(end-7,end-3) + by3_bm;
M(end-6,end-5) = M(end-6,end-5) + bz1_bm;
M(end-6,end-4) = M(end-6,end-4) + bz2_bm;
M(end-6,end-3) = M(end-6,end-3) + bz3_bm;


% Maxwell's equations
idx_Jx_bm = getIndexOfBlkMatrix_maxwell(S,Ns,J,1);
M(end-5,end-1) = M(end-5,end-1) + c2*kz;
M(end-5,end-8) = M(end-5,end-8) - 1i/epsilon0;
M(end-5,idx_Jx_bm(1:end-1)) = M(end-5,idx_Jx_bm(1:end-1)) - 1i/epsilon0;
idx_Jy_bm = getIndexOfBlkMatrix_maxwell(S,Ns,J,2);
M(end-4,end-2) = M(end-4,end-2) - c2*kz;
M(end-4,end) = M(end-4,end) + c2*kx;
M(end-4,end-7) = M(end-4,end-7) - 1i/epsilon0; 
M(end-4,idx_Jy_bm(1:end-1)) = M(end-4,idx_Jy_bm(1:end-1)) - 1i/epsilon0;
idx_Jz_bm = getIndexOfBlkMatrix_maxwell(S,Ns,J,3);
M(end-3,end-1) = M(end-3,end-1) - c2*kx;
M(end-3,end-6) = M(end-3,end-6) - 1i/epsilon0;
M(end-3,idx_Jz_bm(1:end-1)) = M(end-3,idx_Jz_bm(1:end-1)) - 1i/epsilon0;
M(end-2,end-4) = M(end-2,end-4) - kz;
M(end-1,end-5) = M(end-1,end-5) + kz; 
M(end-1,end-3) = M(end-1,end-3) - kx;
M(end,end-4) = M(end,end-4) + kx;

% Solver
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