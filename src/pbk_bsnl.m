function bsnl = pbk_bsnl(s,n,l,kappasz,vtsz,wcs,kz)
% @Description: Calculate the coefficients of bsnl for the oblique 
% plasma waves with a loss-cone product-bi-kappa distribution (for PBK of type 2).
% @Filename: pbk_logbsnl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-11-09

log_csl = @(s,l) gammaln(kappasz(s)+1) + gammaln(2*kappasz(s)-l+2) ...
    -gammaln(kappasz(s)-l+2) - gammaln(2*kappasz(s)+1) ...
    +(l-1)*log(2i*sqrt(kappasz(s))*kz*vtsz(s));

bsnl = -1.0*n*wcs(s) * exp(log_csl(s,l));

end