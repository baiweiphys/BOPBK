function by32snl = pbk_by32snl(s,n,l,theta,csn,b4snl,wps,wcs,norm_omega)
% PBK_BY32SNL Calculate the coefficient by32snl.
% @Description: Calculate the coefficients of by32snl for the y-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_by32snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps and wcs by NORMALIZED frequency: 
%   wps = wps/norm_omega
%   wcs = wcs/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;
norm_wcs = wcs/norm_omega;

by32snl = -1.0*epsilon0*tan(theta)*norm_wps(s).^2 ...
    *b4snl(s,n,l)*(csn(s,n)-n*norm_wcs(s))/norm_wcs(s);

end
