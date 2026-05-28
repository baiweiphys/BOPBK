function bz31snl = pbk_bz31snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs,norm_omega)
% PBK_BZ31SNL Calculate the coefficient bz31snl.
% @Description: Calculate the coefficients of bz31snl for the z-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bz31snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps and wcs by NORMALIZED frequency: 
%   wps = wps/norm_omega
%   wcs = wcs/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

bz31snl = -1i*epsilon0*tan(theta)^2*norm_wps(s).^2 ...
          *(b1snl(s,n,l)*(csn(s,n)-n*norm_wcs(s))^2 + ...
          2*b2snl(s,n,l)*(csn(s,n)-n*norm_wcs(s)))/norm_wcs(s)^2;

end
