function bz3snj = mk_bz3snj(s,n,jj,theta,b12snj,csnj,wps,wcs,norm_omega)
% MK_BZ3SNJ Calculate the coefficients of bz3snj.
% @Description: Calculate the coefficients of bz3snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_bz3snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps and wcs by NORMALIZED frequency: 
%                        wps = wps/norm_omega, wcs = wcs/norm_omega
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

coef = norm_wps(s)^2*b12snj(s,n,jj)*(csnj(s,n,jj)/norm_wcs(s)^2 - 2*n/norm_wcs(s) + n^2/csnj(s,n,jj));
bz3snj = -1i*epsilon0*tan(theta)^2*coef;

end
