function b34snj = mk1_b34snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,norm_omega,EPS0)
% MK1_B34SNJ Calculate the coefficients of b34snj.
% @Description: Calculate the coefficients of b34snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK, PBK type1) distribution.
% @Filename: mk1_b34snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale b34snj by NORMALIZED frequency: b34snj = b34snj/norm_omega.
%   Bai Wei, 2026.04.02: Coefficient is now conditionally switched:
%                        MK coefficient (default) OR 
%                        BM coefficient if (kappaxs > kappaxs_th)


%% FuncS_pbk = @(Jnum,dJnum,num,den) funcS_pbk1(s,n,kappas_perp,lambdas,sgms,Jnum,dJnum,num,den,pbk_type,EPS0)
% S1_pbk = FuncS_pbk1(2,0,1,2)
% S2_pbk = FuncS_pbk1(2,0,1,1)
% S3_pbk = FuncS_pbk11(1,1,2,2)
% S4_pbk = FuncS_pbk11(1,1,2,1)
% S5_pbk = FuncS_pbk1(0,2,3,2)
% S6_pbk = FuncS_pbk1(0,2,3,1)
% S7_pbk = FuncS_pbk1(2,0,-1,1)
% S8_pbk = FuncS_pbk1(1,1,0,1)
% S9_pbk = FuncS_pbk1(0,2,1,1)

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
    S3_pbk = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,1,1,2,2,EPS0);
    S4_pbk = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,1,1,2,1,EPS0);
    S8_pbk = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,1,1,0,1,EPS0);
    S38_pbk = @(s,n) (kappasx(s)+sgms(s)+1)/kappasx(s)*S3_pbk(s,n)...
        -2*sgms(s)*lambdas(s)*S8_pbk(s,n);
    %
    b34snj_mk = bj(jj)*n*wcs(s)*S38_pbk(s,n) ...
            + bj(jj)*cj(jj)*kz*vtsz(s)*S4_pbk(s,n)*(1-1/kappasx(s))*Tsx(s)/Tsz(s)/(1+sgms(s));
    b34snj = b34snj_mk/norm_omega;
else
    % bi-Maxwellian (BM)
    if sgms(s)==0
        b34snj_tmp = -1*bsnj(s,n,jj)*(In(s,n)-dIn(s,n));
    else
        S3_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,1,1,2,EPS0); %S3_mx=S4_mx 
        S8_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,1,1,0,EPS0);
        S38_max = @(s,n) S3_max(s,n) - 2*sgms(s)*lambdas(s)*S8_max(s,n);

        b34snj_tmp = bj(jj)*n*wcs(s)*S38_max(s,n) ...
            + bj(jj)*cj(jj)*kz*vtsz(s)*S3_max(s,n)*Tsx(s)/Tsz(s)/(1+sgms(s));
    end
    b34snj = b34snj_tmp/norm_omega;
end
