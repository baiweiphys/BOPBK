function bx22snl = pbk_bx22snl(s,n,l,b4snl,wps,norm_omega)
% PBK_BX22SNL Calculate the coefficient bx22snl.
% @Description: Calculate the coefficients of bx22snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx22snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

bx22snl = epsilon0*norm_wps(s).^2*n*b4snl(s,n,l);

end
