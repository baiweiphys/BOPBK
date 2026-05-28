function idx = SNLJmap2oneDim_pbk2(Ns_pbk,kappasz_pbk,s,idx_n,l,j)
% SNLJMAP2ONEDIM_PBK2 Calculate the 1-dimensional PBK array index. 
% @Description: Calculate the 1-dimensional PBK array index corresponding 
% to the given parameters (N_pbk2,kappasz_pbk2,s,idx_n,l,j)
% @Filename: SNLJmap2oneDim_pbk2.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2025-02-06
% @Modification: 
%   Bai Wei, 2026.04.04: update SNLJmap2oneDim for pbk type2.


 if idx_n < 1 || idx_n > 2 * Ns_pbk(s) + 1
     error('n must be in the range 1 to 2 * N(s) + 1');
 end
 if l < 1 || l > kappasz_pbk(s)
     % PBK type 2
     error('l must be in the range 1 to K(s)');
 end   
 if j < 1 || j > l + 1
     error('j must be in the range 1 to l + 1');
 end


  % 1. Count the elements at the species s-1.
  soff_pbk2 = 0;
  for ii = 1:(s-1)
      tmp = (2*Ns_pbk(ii) + 1)*kappasz_pbk(ii)*(kappasz_pbk(ii)+3)/2;
      soff_pbk2  = soff_pbk2 + tmp;
  end

  % 2. Count the elements at the index n-1.  
  noff_pbk2 = (idx_n-1)*(kappasz_pbk(s)+3)*kappasz_pbk(s)/2;

  % 3. Calculate the offset for the index l-1.
  loff_pbk2 = (l + 2)*(l - 1)/2;

  idx = soff_pbk2 + noff_pbk2  + loff_pbk2 + j;
  % fprintf('1-dimensional array index: %d\n', index);
end
