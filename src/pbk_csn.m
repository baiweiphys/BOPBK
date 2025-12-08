function csn = pbk_csn(s,n,kz,kappasz_pbk,vtsz_pbk,wcs,us0)
% @Description: Calculate the coefficients of csn(s) for the oblique plasma 
% waves with a loss-cone PBK distribution.
% @Filename: csn.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-06

csn = n*wcs(s) + kz*us0(s) -1i*sqrt(kappasz_pbk(s))*kz*vtsz_pbk(s);
