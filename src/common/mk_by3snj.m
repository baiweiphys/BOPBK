function by3snj = mk_by3snj(s,n,jj,theta,b34snj,csnj,wps,wcs,norm_omega)
% MK_BY3SNJ Calculate the coefficients of by3snj.
% @Description: Calculate the coefficients of by3snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_by3snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps and wcs by NORMALIZED frequency: 
%                        wps = wps/norm_omega, wcs = wcs/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;


by3snj = -1*epsilon0*tan(theta)*norm_wps(s)^2*(1-n*norm_wcs(s)/csnj(s,n,jj))*b34snj(s,n,jj)/norm_wcs(s);

end
