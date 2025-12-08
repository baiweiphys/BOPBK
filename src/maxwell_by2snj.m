function by2snj = maxwell_by2snj(s,n,jj,b56snj,csnj,wps)
% @Description: Calculate the coefficients of by2snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_by2snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

by2snj = -1i*epsilon0*wps(s)^2*b56snj(s,n,jj)/csnj(s,n,jj);

end