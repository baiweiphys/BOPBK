function b34snj = maxwell_b34snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0)
% @Description: Calculate the coefficients of b34snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_b34snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-06

% FuncS_max = @(Jnum,dJnum,num) funcS_max(s,n,lambdas,sgms,Jnum,dJnum,num,EPS0)
% S1_max = S2_max = FuncS_max(2,0,1)
% S3_max = S4_max = FuncS_max(1,1,2)
% S5_max = S6_max = FuncS_max(0,2,3)
% S7_max = FuncS_max(2,0,-1)
% S8_max = FuncS_max(1,1,0)
% S9_max = FuncS_max(0,2,1)


% % J-pole
[bj,cj] = func_Jpole(J_opt);

bsnj = @(s,n,jj) maxwell_bsnj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,wcs);

scale = 1;
In = @(s,n) besseli(n,lambdas(s),scale);
% fixed bug: besseli(n+1,lambdas(s),scale), 2024.04.28
dIn = @(s,n) 0.5*besseli(n+1,lambdas(s),scale) + 0.5*besseli(n-1,lambdas(s),scale);


if sgms(s)==0
    b34snj = -1*bsnj(s,n,jj)*(In(s,n)-dIn(s,n));
else
    S3_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,1,1,2,EPS0); %S3_mx=S4_mx 
    S8_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,1,1,0,EPS0);
    S38_max = @(s,n) S3_max(s,n) - 2*sgms(s)*lambdas(s)*S8_max(s,n);

    b34snj = bj(jj)*n*wcs(s)*S38_max(s,n) ...
        + bj(jj)*cj(jj)*kz*vtsz(s)*S3_max(s,n)*Tsx(s)/Tsz(s)/(1+sgms(s));
end
