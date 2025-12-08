function bx33snl = pbk_bx33snl(s,n,l,theta,b1snl,wps,wcs)
% @Description: Calculate the coefficients of bx33snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx33snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-2-06

params_with_unit;

bx33snl = -1i*epsilon0*tan(theta)*wps(s).^2*n*b1snl(s,n,l)./wcs(s);

end