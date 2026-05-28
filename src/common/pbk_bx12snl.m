function bx12snl = pbk_bx12snl(s,n,l,b2snl,wps,norm_omega)
% PBK_BX12SNL Calculate the coefficients of bx12snl.
% @Description: Calculate the coefficients of bx12snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx12snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

bx12snl = -1i*epsilon0*norm_wps(s).^2*n^2*b2snl(s,n,l);

end
