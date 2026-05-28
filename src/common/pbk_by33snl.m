function by33snl = pbk_by33snl(s,n,l,theta,b3snl,wps,wcs,norm_omega)
% PBK_BY33SNL Calculate the coefficient by33snl.
% @Description: Calculate the coefficients of by33snl for the y-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_by33snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps and wcs by NORMALIZED frequency: 
%   wps = wps/norm_omega
%   wcs = wcs/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

by33snl = -1.0*epsilon0*tan(theta)*norm_wps(s).^2*b3snl(s,n,l)/norm_wcs(s);

end
