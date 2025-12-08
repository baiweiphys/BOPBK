function b12snj = maxwell_b12snj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,sgms,wcs,lambdas,EPS0)
% @Description: Calculate the coefficients of b12snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_b12snj.m
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


if sgms(s)==0
    b12snj = bsnj(s,n,jj)*In(s,n)/lambdas(s);
else
    S1_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,2,0,1,EPS0); % S1_mx=S2_mx
    S7_max = @(s,n) funcS_maxwell(s,n,lambdas,sgms,2,0,-1,EPS0);
    S17_max = @(s,n) S1_max(s,n) - 2*sgms(s)*lambdas(s)*S7_max(s,n);

    b12snj = bj(jj)*n*wcs(s)*S17_max(s,n) ...
        + bj(jj)*cj(jj)*kz*vtsz(s)*S1_max(s,n)*Tsx(s)/Tsz(s)/(1+sgms(s));
end
