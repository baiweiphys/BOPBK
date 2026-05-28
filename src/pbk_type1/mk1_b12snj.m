function b12snj = mk1_b12snj(s,n,jj,kz,kappasx,kappasx_th,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,norm_omega,EPS0)
% MK1_B12SNJ Calculate the coefficients of b12snj. 
% @Description: Calculate the coefficients of b12snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK, PBK type1) distribution.
% @Filename: mk1_b12snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale b12snj by NORMALIZED frequency: b12snj = b12snj/norm_omega.
%   Bai Wei, 2026.04.02: Coefficient is now conditionally switched:
%                        MK coefficient (default) OR 
%                        BM coefficient if (kappaxs > kappaxs_th)


%% FuncS_pbk = @(Jnum,dJnum,num,den) funcS_pbk1(s,n,kappas_perp,lambdas,sgms,Jnum,dJnum,num,den,pbk_type,EPS0)
% S1_pbk = FuncS_pbk1(2,0,1,2)
% S2_pbk = FuncS_pbk1(2,0,1,1)
% S3_pbk = FuncS_pbk1(1,1,2,2)
% S4_pbk = FuncS_pbk1(1,1,2,1)
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

%%
if kappasx(s)<kappasx_th(s)
    % Maxwellian-kappa (MK)
    S1_pbk = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,2,0,1,2,EPS0);
    S2_pbk = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,2,0,1,1,EPS0);
    S7_pbk = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,2,0,-1,1,EPS0);
    S17_pbk = @(s,n) (kappasx(s)+sgms(s)+1)/kappasx(s)*S1_pbk(s,n) ...
        -2*sgms(s)*lambdas(s)*S7_pbk(s,n);
    %
    b12snj_mk = bj(jj)*n*wcs(s)*S17_pbk(s,n) ...
            + bj(jj)*cj(jj)*kz*vtsz(s)*S2_pbk(s,n)*(1-1/kappasx(s))*Tsx(s)/Tsz(s)/(1+sgms(s));
    b12snj = b12snj_mk/norm_omega;
else
    % bi-Maxwellian (BM)
    if sgms(s)==0
        b12snj_bm = bsnj(s,n,jj)*In(s,n)/lambdas(s);
    else
        S1_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,2,0,1,EPS0); % S1_mx=S2_mx
        S7_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,2,0,-1,EPS0);
        S17_max = @(s,n) S1_max(s,n) - 2*sgms(s)*lambdas(s)*S7_max(s,n);

        b12snj_bm = bj(jj)*n*wcs(s)*S17_max(s,n) ...
            + bj(jj)*cj(jj)*kz*vtsz(s)*S1_max(s,n)*Tsx(s)/Tsz(s)/(1+sgms(s));
    end
    b12snj = b12snj_bm/norm_omega;
end

