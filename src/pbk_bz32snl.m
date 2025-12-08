function bz32snl = pbk_bz32snl(s,n,l,theta,csn,b2snl,wps,wcs)
% @Description: Calculate the coefficients of bz32snl for the z-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bz32snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bz32snl = -1i*epsilon0*tan(theta)^2*wps(s).^2*b2snl(s,n,l)*(csn(s,n)-n*wcs(s))^2/wcs(s)^2;

end