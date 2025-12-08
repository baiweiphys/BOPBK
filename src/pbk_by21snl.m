function by21snl = pbk_by21snl(s,n,l,b5snl,wps)
% @Description: Calculate the coefficients of by21snl for the y-component 
% of oblique plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_by12snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-1-19

params_with_unit;

by21snl = -1i*epsilon0*wps(s).^2*b5snl(s,n,l);

end