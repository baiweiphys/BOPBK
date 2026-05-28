function bx21snl = pbk_bx21snl(s,n,l,b3snl,wps,norm_omega)
% PBK_BX21SNL Calculate the coefficient bx21snl.
% @Description: Calculate the coefficients of bx21snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx21snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

bx21snl = epsilon0*norm_wps(s).^2*n*b3snl(s,n,l);

end
