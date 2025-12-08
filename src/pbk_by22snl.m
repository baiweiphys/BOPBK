function by22snl = pbk_by22snl(s,n,l,b6snl,wps)
% @Description: Calculate the coefficients of by22snl for the y-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_by22snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

by22snl = -1i*epsilon0*wps(s).^2*b6snl(s,n,l);

end