function idx = SNJmap2oneDim(Ns_bm,J,s,idx_n,j)
% @Description: Calculate the 1-dimensional bi-Maxwellian array index corresponding 
% to the given parameters (N_max,J,s,idx_n,j)
% @Filename: SNJmap2oneDim.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2025-02-06
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-06

if (idx_n < 1 || idx_n > 2*Ns_bm(s) + 1)
    error('n must be in the range 1 to 2 * N(s) + 1');
end
if (j < 1 || j > J)
    error('j must be in the range 1 to J');
end

idx = 0;
for i = 1:(s-1)
    totalElements_i = (2 * Ns_bm(i) + 1) * J;
    idx = idx + totalElements_i;
end
for i = 1:(idx_n-1)
    idx = idx + J;
end
idx = idx + j;
% fprintf('1-dimensional array index: %d\n', index);

end
