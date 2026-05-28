function bx11snl = pbk_bx11snl(s,n,l,b1snl,wps,norm_omega)
% PBK_BX11SNL Calculate the coefficient bx11snl.
% @Description: Calculate the coefficients of bx11snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx11snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

bx11snl = -1i*epsilon0*norm_wps(s).^2*n^2*b1snl(s,n,l);

end
