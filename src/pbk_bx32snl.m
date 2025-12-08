function bx32snl = pbk_bx32snl(s,n,l,theta,csn,b2snl,wps,wcs)
% @Description: Calculate the coefficients of bx32snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx32snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bx32snl = -1i*epsilon0*tan(theta)*wps(s).^2*n*b2snl(s,n,l)*(csn(s,n)-n*wcs(s))/wcs(s);

end