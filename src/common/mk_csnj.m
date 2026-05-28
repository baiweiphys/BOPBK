function csnj = mk_csnj(s,n,jj,kz,J_opt,vtsz,wcs,us0,norm_omega)
% MK_CSNJ Calculate the coefficient csnj.
% @Description: Calculate the coefficient csnj for the oblique plasma 
% wave model with a Maxwellian-kappa (MK) distribution.
% @Filename: mk_csnj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale csnj by NORMALIZED frequency: csnj = csnj/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

[~,cj] = func_Jpole(J_opt);
csnj_tmp = n*wcs(s) + kz*us0(s) + cj(jj)*kz*vtsz(s);

csnj = csnj_tmp/norm_omega;

end
