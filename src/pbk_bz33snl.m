function bz33snl = pbk_bz33snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs)
% @Description: Calculate the coefficients of bz33snl for the z-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bz33snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bz33snl = -1i*epsilon0*tan(theta)^2*wps(s).^2*(2*b1snl(s,n,l)*(csn(s,n)-n*wcs(s)) + b2snl(s,n,l))/wcs(s)^2;

end