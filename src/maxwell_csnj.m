function csnj = maxwell_csnj(s,n,jj,kz,J_opt,vtsz,wcs,us0)
% @Description: Calculate the coefficients for the oblique plasma 
% wave model with a bi-Maxwellian distribution.
% @Filename: maxwell_csnj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

[~,cj] = func_Jpole(J_opt);
csnj = n*wcs(s) + kz*us0(s) + cj(jj)*kz*vtsz(s);

end
