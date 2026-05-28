function by1snj = mk_by1snj(s,n,j,b34snj,csnj,wps,norm_omega)
% MK_BY1SNJ Calculate the coefficients of by1snj.
% @Description: Calculate the coefficients of by1snj for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_by1snj.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps by NORMALIZED frequency: wps = wps/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

params_with_unit;
norm_wps = wps/norm_omega;

by1snj = -1*epsilon0*norm_wps(s)^2*n*b34snj(s,n,j)/csnj(s,n,j);

end
