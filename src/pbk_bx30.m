function bx30 = pbk_bx30(theta,wps)
% @Description: Calculate the coefficients of bx30 for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx30.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bx30 = -1i*epsilon0*tan(theta)*sum(wps.^2);

end