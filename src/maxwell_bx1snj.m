function bx1snj = maxwell_bx1snj(s,n,jj,b12snj,csnj,wps)
% @Description: Calculate the coefficients of bx1snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_bx1snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

bx1snj = -1i*epsilon0*wps(s)^2*n^2*b12snj(s,n,jj)/csnj(s,n,jj);

end