function bx31snl = pbk_bx31snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs,norm_omega)
% PBK_BX31SNL Calculate the coefficient bx31snl.
% @Description: Calculate the coefficients of bx31snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx31snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps and wcs by NORMALIZED frequency: 
%   wps = wps/norm_omega
%   wcs = wcs/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

bx31snl = -1i*epsilon0*tan(theta)*norm_wps(s).^2 ...
    *n*(b1snl(s,n,l)*(csn(s,n)-n*norm_wcs(s))+b2snl(s,n,l))/norm_wcs(s);

end
