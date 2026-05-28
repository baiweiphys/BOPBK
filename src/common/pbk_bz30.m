function bz30 = pbk_bz30(theta,wps,norm_omega)
% PBK_BZ30 Calculate the coefficients bz30.
% @Description: Calculate the coefficients of bz30 for the z-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bz30.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

bz30 = 1i*epsilon0*sum(norm_wps.^2)*tan(theta).^2;
end
