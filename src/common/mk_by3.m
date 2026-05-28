function by3 = mk_by3(S_mk,Ns_mk,J,theta,b34snj,csnj,wps,norm_omega)
% MK_BY3 Calculate the coefficients of by3.
% @Description: Calculate the coefficients of by3 for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_by3.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps by NORMALIZED frequency: wps = wps/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

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

by3 = -1*epsilon0*tan(theta)*coef;

end
