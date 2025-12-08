function bx12snl = pbk_bx12snl(s,n,l,b2snl,wps)
% @Description: Calculate the coefficients of bx12snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx12snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;
bx12snl = -1i*epsilon0*wps(s).^2*n^2*b2snl(s,n,l);

end