function bz1snj = maxwell_bz1snj(s,n,jj,theta,b12snj,csnj,wps,wcs)
% @Description: Calculate the coefficients of bz1snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_bz1snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

coef = wps(s)^2*n*(1-n*wcs(s)/csnj(s,n,jj))*b12snj(s,n,jj)/wcs(s);
bz1snj = -1i*epsilon0*tan(theta)*coef;

end