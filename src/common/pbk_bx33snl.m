function bx33snl = pbk_bx33snl(s,n,l,theta,b1snl,wps,wcs,norm_omega)
% PBK_BX33SNL Calculate the coefficient bx33snl.
% @Description: Calculate the coefficients of bx33snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx33snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps and wcs by NORMALIZED frequency: 
%   wps = wps/norm_omega
%   wcs = wcs/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

bx33snl = -1i*epsilon0*tan(theta)*norm_wps(s).^2*n*b1snl(s,n,l)./norm_wcs(s);

end
