function by22snl = pbk_by22snl(s,n,l,b6snl,wps,norm_omega)
% PBK_BY22SNL Calculate the coefficient by22snl.
% @Description: Calculate the coefficients of by22snl for the y-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_by22snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

by22snl = -1i*epsilon0*norm_wps(s).^2*b6snl(s,n,l);

end
