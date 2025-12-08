function bx21snl = pbk_bx21snl(s,n,l,b3snl,wps)
% @Description: Calculate the coefficients of bx21snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx21snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bx21snl = epsilon0*wps(s).^2*n*b3snl(s,n,l);

end