function bx10 = pbk_bx10(wps,kappasx,sgms)
% @Description: Calculate the coefficients of bx10 for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx10.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-01-19

params_with_unit;

% BUGFIX: Removed redundant coefs (baiwei 2025-06-16)
bx10 = 1i*epsilon0*sum(wps.^2);

end