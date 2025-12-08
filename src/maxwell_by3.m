function by3 = maxwell_by3(S_bm,Ns_bm,J,theta,b34snj,csnj,wps)
% @Description: Calculate the coefficients of by3 for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian distribution.
% @Filename: maxwell_by3.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13

params_with_unit;

coef = 0;
for s=1:S_bm
    Nvec = -Ns_bm(s):Ns_bm(s);
    for idx_n=1:(2*Ns_bm(s)+1)
        n = Nvec(idx_n);
        for jj=1:J
            coef = coef + wps(s)^2*n*b34snj(s,n,jj)/csnj(s,n,jj);
        end
    end
end

by3 = -1*epsilon0*tan(theta)*coef;

end