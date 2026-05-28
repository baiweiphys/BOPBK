function idx = SNLJmap2oneDim_pbk1(Ns_pbk,kappasz_pbk,s,ind_n,l,j)
% SNLJMAP2ONEDIM_PBK1 Calculate the 1-dimensional PBK (type 1) array index. 
% @Description: Calculate the 1-dimensional PBK (type 1) array index corresponding 
% to the given parameters (N_pbk,kappasz_pbk,s,ind_n,l,j)
% @Filename: SNLJmap2oneDim_pbk1.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2025-02-06
% @Modification: 
%   Bai Wei, 2026.04.02: only updated to pbk1.

 if ind_n < 1 || ind_n > 2 * Ns_pbk(s) + 1
     error('n must be in the range 1 to 2 * N(s) + 1');
 end
 if l < 1 || l > kappasz_pbk(s)+1
     % PBK type 1
     error('l must be in the range 1 to K(s)+1');
 end   
 if j < 1 || j > l + 1
     error('j must be in the range 1 to l + 1');
 end


  % 1. Count the elements at the species s-1.
  soff_pbk1 = 0;
  for ii = 1:(s-1)
      tmp = (2*Ns_pbk(ii) + 1)*(kappasz_pbk(ii)+1)*(kappasz_pbk(ii)+4)/2;
      soff_pbk1  = soff_pbk1 + tmp;
  end

  % 2. Count the elements at the index n-1.
  noff_pbk1 = (ind_n-1)*(kappasz_pbk(s)+4)*(kappasz_pbk(s)+1)/2;

  % 3. Calculate the offset for the index l-1.
  loff_pbk1 = (l + 2)*(l - 1)/2;

   idx = soff_pbk1 + noff_pbk1  + loff_pbk1 + j;

    % fprintf('1-dimensional array index: %d\n', index);
end
