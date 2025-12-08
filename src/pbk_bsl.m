function bsl = pbk_bsl(s,l,kappasz,vtsz,vtsx,kz)
% @Description: Calculate the coefficients of bsl for the oblique 
% plasma waves with a product-bi-kappa distribution (for PBK of type 2).
% @Filename: pbk_logbsl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @LastEditors: Bai Wei
% @LastEditTime: 2025-11-09

log_csl = @(s,l) gammaln(kappasz(s)+1) + gammaln(2*kappasz(s)-l+2) ...
    -gammaln(kappasz(s)-l+2) - gammaln(2*kappasz(s)+1) ...
    +(l-1)*log(2i*sqrt(kappasz(s))*kz*vtsz(s));

bsl = -0.5*l*kz^2.*vtsx(s).^2 * exp(log_csl(s,l));

end
