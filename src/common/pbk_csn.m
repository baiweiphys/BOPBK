function csn = pbk_csn(s,n,kz,kappasz,vtsz,wcs,us0,norm_omega)
% PBK_CSN Calculate the coefficient csn.
% @Description: Calculate the coefficients of csn for the oblique plasma 
% waves with a loss-cone PBK distribution.
% @Filename: pbk_csn.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale wcs, kz*us0, and kz*vtsz(s) 
%   by NORMALIZED frequency:
%   wcs = wcs/norm_omega
%   us0 = kz*us0/norm_omega
%   kz*vtsz = kz*vtsz/norm_omega;

norm_wcs = wcs/norm_omega;
norm_us0 = kz*us0/norm_omega;
norm_vtsz = kz*vtsz/norm_omega;

csn = n*norm_wcs(s) + norm_us0(s) -1i*sqrt(kappasz(s))*norm_vtsz(s);
