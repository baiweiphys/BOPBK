function bx10 = pbk_bx10(wps,kappasx,sgms,norm_omega)
% PBK_BX10 Calculate the coefficient bx10.
% @Description: Calculate the coefficients of bx10 for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx10.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

% BUGFIX: Removed redundant coefs (baiwei 2025-06-16)
bx10 = 1i*epsilon0*sum(norm_wps.^2);

end
