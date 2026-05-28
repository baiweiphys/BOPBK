function by2snj = mk_by2snj(s,n,jj,b56snj,csnj,wps,norm_omega)
% MK_BY2SNJ Calculate the coefficients of by2snj.
% @Description: Calculate the coefficients of by2snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_by2snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps by NORMALIZED frequency: wps = wps/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

params_with_unit;
norm_wps = wps/norm_omega;

by2snj = -1i*epsilon0*norm_wps(s)^2*b56snj(s,n,jj)/csnj(s,n,jj);

end
