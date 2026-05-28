function bsl = pbk1_bsl(s,l,kappasz,vtsz,vtsx,kz,norm_omega)
% PBK1_BSL Calculate the coefficient bsl.
% @Description: Calculate the coefficients of bsl for the oblique 
% plasma waves with a product-bi-kappa distribution (PBK type 1/2 model).
% @Filename: pbk1_bsl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-24
% @Modification: 
%   Bai Wei, 2026.04.02: Scale kz*vtsz(s) and kz*vtsx(s) by NORMALIZED frequency: 
%   kz*vtsz = kz*vtsz/norm_omega;
%   kx*vtsx = kz*vtsx/norm_omega

norm_vtsz = kz*vtsz/norm_omega;
norm_vtsx = kz*vtsx/norm_omega;

% PBK type 1
log_csl = @(s,l) gammaln(kappasz(s)+1) + gammaln(2*kappasz(s)-l+2) ...
    -gammaln(kappasz(s)-l+2) - gammaln(2*kappasz(s)+1) ...
    +(l-1)*log(2i*sqrt(kappasz(s))*norm_vtsz(s));

bsl = -0.5*l*norm_vtsx(s).^2 * exp(log_csl(s,l));

end
