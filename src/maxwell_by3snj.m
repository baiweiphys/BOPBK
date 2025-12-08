function by3snj = maxwell_by3snj(s,n,jj,theta,b34snj,csnj,wps,wcs)
% @Description: Calculate the coefficients of by3snj for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_by3snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

by3snj = -1*epsilon0*tan(theta)*wps(s)^2*(1-n*wcs(s)/csnj(s,n,jj))*b34snj(s,n,jj)/wcs(s);

end
