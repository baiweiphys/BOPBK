function bz2snj = mk_bz2snj(s,n,jj,theta,b34snj,csnj,wps,wcs,norm_omega)
% MK_BZ2SNJ Calculate the coefficients of bz2snj.
% @Description: Calculate the coefficients of bz2snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_bz2snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps and wcs by NORMALIZED frequency: 
%                        wps = wps/norm_omega, wcs = wcs/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

coef = norm_wps(s)^2*(1-n*norm_wcs(s)/csnj(s,n,jj))*b34snj(s,n,jj)/norm_wcs(s);
bz2snj = epsilon0*tan(theta)*coef;

end
