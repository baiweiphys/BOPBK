function b56snj = mk2_b56snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,norm_omega,EPS0)
% MK_B56SNJ Calculate the coefficients of b56snj.
% @Description: Calculate the coefficients of b56snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK, PBK type 2) distribution.
% @Filename: mk2_b56snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Modification: 
%   Bai Wei, 2026.04.01: Scale b56snj by NORMALIZED frequency: b56snj = b56snj/norm_omega.
%   Bai Wei, 2026.04.02: Coefficient is now conditionally switched:
%                        MK coefficient (default) OR 
%                        BM coefficient if (kappaxs > kappaxs_th)
%   Bai Wei, 2026.04.04: update S3, S4, S8 and S38 for pbk type2.


%% FuncS_pbk = @(Jnum,dJnum,num,den) funcS_pbk2(s,n,kappas_perp,lambdas,sgms,Jnum,dJnum,num,den,EPS0)
% S1_pbk = FuncS_pbk2(2,0,1,1)
% S2_pbk = FuncS_pbk2(2,0,1,0)
% S3_pbk = FuncS_pbk2(1,1,2,1)
% S4_pbk = FuncS_pbk2(1,1,2,0)
% S5_pbk = FuncS_pbk2(0,2,3,1)
% S6_pbk = FuncS_pbk2(0,2,3,0)
% S7_pbk = FuncS_pbk2(2,0,-1,0)
% S8_pbk = FuncS_pbk2(1,1,0,0)
% S9_pbk = FuncS_pbk2(0,2,1,0)

%% FuncS_max = @(Jnum,dJnum,num) funcS_max(s,n,lambdas,sgms,Jnum,dJnum,num,EPS0)
% S1_max = S2_max = FuncS_max(2,0,1)
% S3_max = S4_max = FuncS_max(1,1,2)
% S5_max = S6_max = FuncS_max(0,2,3)
% S7_max = FuncS_max(2,0,-1)
% S8_max = FuncS_max(1,1,0)
% S9_max = FuncS_max(0,2,1)

%% J-pole
[bj,cj] = func_Jpole(J_opt);

%% for bi-Maxwellian (BM)
bsnj = @(s,n,jj) maxwell_bsnj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,wcs);

scale = 1;
In = @(s,n) besseli(n,lambdas(s),scale);
% fixed bug: besseli(n+1,lambdas(s),scale), 2024.04.28
dIn = @(s,n) 0.5*besseli(n+1,lambdas(s),scale) + 0.5*besseli(n-1,lambdas(s),scale);

%%
if kappasx(s)<kappasx_th(s)
    % Maxwellian-kappa (MK)
    S5_pbk = @(s,n) funcS_pbk2(s,n,kappasx,lambdas,sgms,0,2,3,1,EPS0);
    S6_pbk = @(s,n) funcS_pbk2(s,n,kappasx,lambdas,sgms,0,2,3,0,EPS0);
    S9_pbk = @(s,n) funcS_pbk2(s,n,kappasx,lambdas,sgms,0,2,1,0,EPS0);
    S59_pbk = @(s,n) (kappasx(s)+sgms(s))/kappasx(s)*S5_pbk(s,n) ...
    -2*sgms(s)*lambdas(s)*S9_pbk(s,n);
    %
    b56snj_mk = bj(jj)*n*wcs(s)*S59_pbk(s,n) ...
            + bj(jj)*cj(jj)*kz*vtsz(s)*S6_pbk(s,n)*(1-1/kappasx(s))*Tsx(s)/Tsz(s)/(1+sgms(s));
    b56snj = b56snj_mk/norm_omega;
else
    % bi-Maxwellian (BM)
    if sgms(s)==0
        b56snj_tmp = bsnj(s,n,jj)*(n^2*In(s,n) + 2*lambdas(s)^2*(In(s,n)-dIn(s,n)))/lambdas(s);
    else
        S5_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,0,2,3,EPS0); % S5_mx=S6_mx
        S9_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,0,2,1,EPS0);
        S59_max = @(s,n) S5_max(s,n) - 2*sgms(s)*lambdas(s)*S9_max(s,n);

        b56snj_tmp = bj(jj)*n*wcs(s)*S59_max(s,n) ...
            + bj(jj)*cj(jj)*kz*vtsz(s)*S5_max(s,n)*Tsx(s)/Tsz(s)/(1+sgms(s));
    end
    b56snj = b56snj_tmp/norm_omega;
end
