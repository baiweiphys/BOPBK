function bx32snl = pbk_bx32snl(s,n,l,theta,csn,b2snl,wps,wcs,norm_omega)
% PBK_BX32SNL Calculate the coefficient bx32snl.
% @Description: Calculate the coefficients of bx32snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx32snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps and wcs by NORMALIZED frequency: 
%   wps = wps/norm_omega
%   wcs = wcs/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

bx32snl = -1i*epsilon0*tan(theta)*norm_wps(s).^2*...
    n*b2snl(s,n,l)*(csn(s,n)-n*norm_wcs(s))/norm_wcs(s);

end
