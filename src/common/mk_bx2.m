function bx2 = mk_bx2(S_mk,Ns_mk,J,b34snj,csnj,wps,norm_omega)
% MK_BX2 Calculate the coefficients of bx2.
% @Description: Calculate the coefficients of bx2 for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_bx2.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps by NORMALIZED frequency: wps = wps/norm_omega.
%   Bai Wei, 2026.04.02: Coefficient is now conditionally switched:
%                        MK coefficient (default) OR 
%                        BM coefficient if (kappaxs > kappaxs_th)


params_with_unit;
norm_wps = wps/norm_omega;

coef = 0;
for s=1:S_mk
    Nvec = -Ns_mk(s):Ns_mk(s);
    for ind_n=1:(2*Ns_mk(s)+1)
        n = Nvec(ind_n);
        for jj=1:J
            coef = coef + norm_wps(s)^2*n*b34snj(s,n,jj)/csnj(s,n,jj);
        end
    end
end
bx2 = -1*epsilon0*coef;

end
