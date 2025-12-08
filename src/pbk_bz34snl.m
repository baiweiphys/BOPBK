function bz34snl = pbk_bz34snl(s,n,l,theta,b1snl,wps,wcs)
% @Description: Calculate the coefficients of bz34snl for the z-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bz34snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bz34snl = -1i*epsilon0*tan(theta)^2*wps(s).^2*b1snl(s,n,l)/wcs(s)^2;

end