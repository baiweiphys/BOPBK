function len = getLen_SNLJmap2oneDim(S_pbk,Ns_pbk,kappasz_pbk)
% @Description: Calculate the len of 1-dimensional PBK array corresponding 
% to the given parameters (S_pbk,N_pbk,kappasz_pbk)
% @Filename: getLen_SNLJmap2oneDim.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2025-02-06
% @LastEditors: Bai Wei
% @LastEditTime: 2025-11-07

len = 0;
for s_val = 1:S_pbk
    N_current = Ns_pbk(s_val);
    kappa_current = kappasz_pbk(s_val);
    for n_val = 1:(2 * N_current + 1)
        for l_val = 1:kappa_current+1
            for j_val = 1:(l_val + 1)
                len = len + 1;
            end
        end
    end
end
% fprintf('length of the 1-dimensional array: %d\n', len);

end
