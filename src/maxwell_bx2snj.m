function bx2snj = maxwell_bx2snj(s,n,jj,b34snj,csnj,wps)
% @Description: Calculate the coefficients of bx2snj for the oblique 
% plasma wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_bx2snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

coef = wps(s)^2*n*b34snj(s,n,jj)/csnj(s,n,jj);
bx2snj = epsilon0*coef;

end