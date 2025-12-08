function len = getLen_SNJmap2oneDim(S_bm,Ns_bm,J)
% @Description: Calculate the len of 1-dimensional bi-Maxwellian array 
% corresponding to the given parameters (S_max,N_max,J)
% @Filename: getLen_SNJmap2oneDim.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2025-02-06
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-06


len = 0;
for s_val = 1:S_bm
    N_current = Ns_bm(s_val);
    for n_val = 1:(2 * N_current + 1)
        for j_val = 1:J
            len = len + 1;
        end
    end
end

% fprintf('length of the 1-dimensional array: %d\n', len);
