function bz31snl = pbk_bz31snl(s,n,l,theta,csn,b1snl,b2snl,wps,wcs)
% @Description: Calculate the coefficients of bz31snl for the z-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bz31snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bz31snl = -1i*epsilon0*tan(theta)^2*wps(s).^2 ...
          *(b1snl(s,n,l)*(csn(s,n)-n*wcs(s))^2 + 2*b2snl(s,n,l)*(csn(s,n)-n*wcs(s)))/wcs(s)^2;

end