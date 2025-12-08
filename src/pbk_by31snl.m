function by31snl = pbk_by31snl(s,n,l,theta,csn,b3snl,b4snl,wps,wcs)
% @Description: Calculate the coefficients of by31snl for the y-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_by31snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

by31snl = -1.0*epsilon0*tan(theta)*wps(s).^2*(b3snl(s,n,l)*(csn(s,n)-n*wcs(s)) + b4snl(s,n,l))/wcs(s);

end