function bz34snl = pbk_bz34snl(s,n,l,theta,b1snl,wps,wcs,norm_omega)
% PBK_BZ34SNL Calculate the coefficient bz34snl.
% @Description: Calculate the coefficients of bz34snl for the z-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bz34snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps and wcs by NORMALIZED frequency: 
%   wps = wps/norm_omega
%   wcs = wcs/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

bz34snl = -1i*epsilon0*tan(theta)^2*norm_wps(s).^2*b1snl(s,n,l)/norm_wcs(s)^2;

end
