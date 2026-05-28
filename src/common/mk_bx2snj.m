function bx2snj = mk_bx2snj(s,n,jj,b34snj,csnj,wps,norm_omega)
% MK_BX2SNJ Calculate the coefficients of bx2snj.
% @Description: Calculate the coefficients of bx2snj for the oblique 
% plasma wave model with a Maxwellian-kappa (MK) distribution.
% @Filename: mk_bx2snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps by NORMALIZED frequency: wps = wps/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

params_with_unit;
norm_wps = wps/norm_omega;

coef = norm_wps(s)^2*n*b34snj(s,n,jj)/csnj(s,n,jj);
bx2snj = epsilon0*coef;

end
