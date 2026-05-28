function by21snl = pbk_by21snl(s,n,l,b5snl,wps,norm_omega)
% PBK_BY21SNL Calculate the coefficient by21snl.
% @Description: Calculate the coefficients of by21snl for the y-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_by21snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

by21snl = -1i*epsilon0*norm_wps(s).^2*b5snl(s,n,l);

end
