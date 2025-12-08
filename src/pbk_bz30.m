function bz30 = pbk_bz30(theta,wps)
% @Description: Calculate the coefficients of bz30 for the z-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bz30.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bz30 = 1i*epsilon0*sum(wps.^2)*tan(theta).^2;
end