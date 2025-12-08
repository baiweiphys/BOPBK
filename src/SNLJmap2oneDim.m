function idx = SNLJmap2oneDim(Ns_pbk,kappasz_pbk,s,idx_n,l,j)
% @Description: Calculate the 1-dimensional PBK array index corresponding 
% to the given parameters (N_pbk,kappasz_pbk,s,idx_n,l,j)
% @Filename: SNLJmap2oneDim.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2025-02-06
% @LastEditors: Bai Wei
% @LastEditTime: 2025-11-07

 if idx_n < 1 || idx_n > 2 * Ns_pbk(s) + 1
     error('n must be in the range 1 to 2 * N(s) + 1');
 end
 if l < 1 || l > kappasz_pbk(s)+1
     error('l must be in the range 1 to K(s)+1');
 end   
 if j < 1 || j > l + 1
     error('j must be in the range 1 to l + 1');
 end

  idx = 0;
  for i = 1:(s-1)
      totalElements_i = 0;
      for n_i = 1:(2 * Ns_pbk(i) + 1)
          for l_i = 1:kappasz_pbk(i)+1
              totalElements_i = totalElements_i + (l_i + 1);
          end
       end
       idx = idx + totalElements_i;
   end

    totalElements_s = 0;
    for n_i = 1:(idx_n-1)
        for l_i = 1:kappasz_pbk(s)+1
            totalElements_s = totalElements_s + (l_i + 1);
        end
    end
    idx = idx + totalElements_s;

    totalElements_n = 0;
    for l_i = 1:(l-1)
        totalElements_n = totalElements_n + (l_i + 1);
    end
    idx = idx + totalElements_n;
    idx = idx + j;
    % fprintf('1-dimensional array index: %d\n', index);
end


