function bz3snj = maxwell_bz3snj(s,n,jj,theta,b12snj,csnj,wps,wcs)
% @Description: Calculate the coefficients of bz3snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_bz3snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;
coef = wps(s)^2*b12snj(s,n,jj)*(csnj(s,n,jj)/wcs(s)^2 - 2*n/wcs(s) + n^2/csnj(s,n,jj));
bz3snj = -1i*epsilon0*tan(theta)^2*coef;

end