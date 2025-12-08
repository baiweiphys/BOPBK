function by1snj = maxwell_by1snj(s,n,j,b34snj,csnj,wps)
% @Description: Calculate the coefficients of by1snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_by1snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

by1snj = -1*epsilon0*wps(s)^2*n*b34snj(s,n,j)/csnj(s,n,j);

end