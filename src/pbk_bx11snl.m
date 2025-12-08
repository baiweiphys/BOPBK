function bx11snl = pbk_bx11snl(s,n,l,b1snl,wps)
% @Description: Calculate the coefficients of bx11snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx11snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;
bx11snl = -1i*epsilon0*wps(s).^2*n^2*b1snl(s,n,l);

end