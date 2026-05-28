function bsnl = pbk2_bsnl(s,n,l,kappasz,vtsz,wcs,kz,norm_omega)
% PBK2_BSNL Calculate the coefficient bsnl.
% @Description: Calculate the coefficients of bsnl for the oblique 
% plasma waves with a loss-cone product-bi-kappa distribution (PBK type 2 model).
% @Filename: pbk2_bsnl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale kz*vtsz(s) and kz*vtsx(s) by NORMALIZED frequency: 
%   kz*vtsz = kz*vtsz/norm_omega;
%   wcs = wcs/norm_omega
%   Bai Wei, 2026.04.04: only update log_csl for pbk type2.

norm_vtsz = kz*vtsz/norm_omega;
norm_wcs = wcs/norm_omega;

% PBK type 2
log_csl = @(s,l) gammaln(kappasz(s)) + gammaln(2*kappasz(s)-l) ...
    - gammaln(2*kappasz(s)-1) - gammaln(kappasz(s)-l+1)  ...
    +(l-1)*log(2i*sqrt(kappasz(s))*norm_vtsz(s));

bsnl = -1.0*n*norm_wcs(s) * exp(log_csl(s,l));

end
