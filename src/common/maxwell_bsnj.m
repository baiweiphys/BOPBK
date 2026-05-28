function bsnj = maxwell_bsnj(s,n,jj,kz,J_opt,vtsz,Tsz,Tsx,wcs)
% MAXWELL_BSNJ Calculate the coefficient bsnl.
% @Description: Calculate the coefficients of bsnl for the oblique 
% plasma wave model with a bi-Maxwellian-kappa (BM) distribution.
% @Filename: maxwell_bsnj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: No normalization needed.

[bj,cj] = func_Jpole(J_opt);
bsnj = bj(jj)*cj(jj)*kz*vtsz(s)*Tsx(s)/Tsz(s) + bj(jj)*n*wcs(s);

end
