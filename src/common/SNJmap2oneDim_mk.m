function idx = SNJmap2oneDim_mk(Ns_mk,J,s,ind_n,j)
% SNJMAP2ONEDIM_MK Calculate the 1-dimensional Maxwellian-kappa array index. 
% @Description: Calculate the 1-dimensional Maxwellian-kappa (MK) array index corresponding 
% to the given parameters (Ns_mk,J,s,ind_n,j)
% @Filename: SNJmap2oneDim_mk.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2025-02-06
% @Modification: 
%   Bai Wei, 2026.04.05: Updated BM to MK.

if (ind_n < 1 || ind_n > 2*Ns_mk(s) + 1)
    error('n must be in the range 1 to 2 * N(s) + 1');
end
if (j < 1 || j > J)
    error('j must be in the range 1 to J');
end

% 1.  Sum of  elements for indeces < s.
soff_mk = 0;
for ii = 1:(s-1)
    tmp = (2 * Ns_mk(ii) + 1) * J;
    soff_mk = soff_mk + tmp;
end

% 2. Count the elements at the idx_n-1
noff_mk = (ind_n-1)*J;

idx = soff_mk  + noff_mk + j;

% fprintf('1-dimensional array index: %d\n', index);
end
