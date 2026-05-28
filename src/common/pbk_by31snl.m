function by31snl = pbk_by31snl(s,n,l,theta,csn,b3snl,b4snl,wps,wcs,norm_omega)
% PBK_BY31SNL Calculate the coefficient by31snl.
% @Description: Calculate the coefficients of by31snl for the y-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_by31snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps and wcs by NORMALIZED frequency: 
%   wps = wps/norm_omega
%   wcs = wcs/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

by31snl = -1.0*epsilon0*tan(theta)*norm_wps(s).^2 ...
    *(b3snl(s,n,l)*(csn(s,n)-n*norm_wcs(s)) + b4snl(s,n,l))/norm_wcs(s);

end
