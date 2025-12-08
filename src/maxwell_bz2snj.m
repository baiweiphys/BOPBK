function bz2snj = maxwell_bz2snj(s,n,jj,theta,b34snj,csnj,wps,wcs)
% @Description: Calculate the coefficients of bz2snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_bz2snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-01-21

params_with_unit;

coef = wps(s)^2*(1-n*wcs(s)/csnj(s,n,jj))*b34snj(s,n,jj)/wcs(s);
bz2snj = epsilon0*tan(theta)*coef;

end