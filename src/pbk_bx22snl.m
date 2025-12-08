function bx22snl = pbk_bx22snl(s,n,l,b4snl,wps)
% @Description: Calculate the coefficients of bx22snl for the x-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_bx22snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

bx22snl = epsilon0*wps(s).^2*n*b4snl(s,n,l);

end