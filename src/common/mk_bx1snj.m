function bx1snj = mk_bx1snj(s,n,jj,b12snj,csnj,wps,norm_omega)
% MK_BX1SNJ Calculate the coefficients of bx1snj.
% @Description: Calculate the coefficients of bx1snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_bx1snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps by NORMALIZED frequency: wps = wps/norm_omega.
%   Bai Wei, 2026.04.02: Coefficient is now conditionally switched:
%                        MK coefficient (default) OR 
%                        BM coefficient if (kappaxs > kappaxs_th)
params_with_unit;

norm_wps = wps/norm_omega;

bx1snj = -1i*epsilon0*norm_wps(s)^2*n^2*b12snj(s,n,jj)/csnj(s,n,jj);

end
