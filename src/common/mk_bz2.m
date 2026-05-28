function bz2 = mk_bz2(S_bm,Ns_bm,J,theta,b34snj,csnj,wps,norm_omega)
% MK_BZ2 Calculate the coefficients of bz2.
% @Description: Calculate the coefficients of bz2 for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_bz2.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps by NORMALIZED frequency: wps = wps/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

params_with_unit;
norm_wps = wps/norm_omega;

coef = 0;
for s=1:S_bm
    Nvec = -Ns_bm(s):Ns_bm(s);
    for ind_n=1:(2*Ns_bm(s)+1)
        n = Nvec(ind_n);
        for jj=1:J
            coef = coef + norm_wps(s)^2*n*b34snj(s,n,jj)/csnj(s,n,jj);
        end
    end
end
bz2 = epsilon0*tan(theta)*coef;

end
