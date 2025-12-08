function bx3 = maxwell_bx3(S_bm,Ns_bm,J,theta,b12snj,csnj,wps)
% @Description: Calculate the coefficients of bx3 for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_bx3.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

coef1 = sum(wps.^2);
coef2 = 0;
%
for s=1:S_bm
    Nvec = -Ns_bm(s):Ns_bm(s);
    for idx_n=1:(2*Ns_bm(s)+1)
        n = Nvec(idx_n);
        for jj=1:J
            coef2 = coef2 + wps(s)^2*n^2*b12snj(s,n,jj)/csnj(s,n,jj);
        end
    end
end

bx3 = -1i*epsilon0*tan(theta)*(coef1 + coef2);

end