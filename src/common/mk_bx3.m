function bx3 = mk_bx3(S_mk,Ns_mk,J,theta,b12snj,csnj,wps,norm_omega)
% MK_BX3 Calculate the coefficients of bx3. 
% @Description: Calculate the coefficients of bx3 for the oblique 
% plasma wave model with a loss-cone Maxwellian-kappa (MK) distribution.
% @Filename: mk_bx3.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.01: Scale wps by NORMALIZED frequency: wps = wps/norm_omega.
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

params_with_unit;
norm_wps = wps/norm_omega;

coef1 = sum(norm_wps.^2);
coef2 = 0;
%
for s=1:S_mk
    Nvec = -Ns_mk(s):Ns_mk(s);
    for ind_n=1:(2*Ns_mk(s)+1)
        n = Nvec(ind_n);
        for jj=1:J
            coef2 = coef2 + norm_wps(s)^2*n^2*b12snj(s,n,jj)/csnj(s,n,jj);
        end
    end
end

bx3 = -1i*epsilon0*tan(theta)*(coef1 + coef2);

end
