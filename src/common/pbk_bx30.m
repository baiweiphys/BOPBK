function bx30 = pbk_bx30(theta,wps,norm_omega)
% PBK_BX30 Calculate the coefficient bx30.
% @Description: Calculate the coefficients of bx30 for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx30.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wps by NORMALIZED frequency: 
%   wps = wps/norm_omega

params_with_unit;
norm_wps = wps/norm_omega;

bx30 = -1i*epsilon0*tan(theta)*sum(norm_wps.^2);

end
