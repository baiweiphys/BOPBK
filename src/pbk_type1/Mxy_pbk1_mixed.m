function M = Mxy_pbk1_mixed(S_pbk,Ns_pbk,Ns_mk,kappasz_pbk,Jpole,csn_pbk, ...
                           bxy11snl,bxy21snl,bxy31snl, ...
                           bxy12snl,bxy22snl,bxy32snl, ...
                           bxy33snl,bx10_by20, ...
                           MatrixNo,ExNo,EyNo,EzNo,ExyNo)
% MXY_PBK1_MIXED Assemble the x- or y-direction PBK block in the mixed solver.
% @Description: To obtain matrix Mxy_pbk_mixed for the mixed loss-cone
% distributions of PBK and bi-Maxwellian plasmas in x- or y-direction.
% @Filename: Mxy_pbk1_mixed.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-10-16
% @Modification: 
%   Bai Wei, 2026.04.02: Updated BM to MK (retains BM features).

% ExyzNo = 5 for Ex
% ExyzNo = 4 for Ey
% ExyzNo = 3 for Ez

% PBK type 1
ind_M_pbk = @(s,ind_n,l,j) SNLJmap2oneDim_pbk1(Ns_pbk,kappasz_pbk,s,ind_n,l,j);

% Step 1
% PBK type 1
len_M_pbk = 0.5*sum((2*Ns_pbk+1).*(kappasz_pbk+4).*(kappasz_pbk+1)) + 1; % PBK submatrix
len_M_mk = Jpole*sum(2*Ns_mk+1) + 1;  % BM submatrix


% Step 2: Create Matrix
% PBK type 1
ind_BlkMatrix = getIndexOfBlkMatrix_mixed_pbk1(Ns_pbk,Ns_mk,kappasz_pbk,Jpole,MatrixNo);
firstInd = ind_BlkMatrix(1)-1;

% for PBK matrix 
len_pbk = 3*len_M_pbk;

% for bi-Maxwellian matrix 
len_bm = 3*len_M_mk;

len_row = len_M_pbk;
len_col = len_pbk + len_bm + 9;
M = zeros(len_row,len_col);


% Step 3: Assemble Matrix
for s=1:S_pbk
    Nvec = -Ns_pbk(s):Ns_pbk(s);
    for ind_n=1:(2*Ns_pbk(s)+1)
        n = Nvec(ind_n);
        for l = 1:kappasz_pbk(s)+1
            for jj = 1:l+1
                snlj = ind_M_pbk(s,ind_n,l,jj);
                M(snlj,firstInd+snlj) = M(snlj,firstInd+snlj) + csn_pbk(s,n); 
                if(jj<l+1) 
                    %snljp1 = index_Mpbk(s,idx_n,l,jj+1);
                    M(snlj,firstInd+snlj+1) = M(snlj,firstInd+snlj+1)+1;
                end
                if (jj==l)
                    if (l<=kappasz_pbk(s))
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bxy33snl(s,n,l+1);  
                        M(snlj,end-ExNo) = M(snlj,end-ExNo) + bxy11snl(s,n,l); 
                        M(snlj,end-EyNo) = M(snlj,end-EyNo) + bxy21snl(s,n,l);
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bxy31snl(s,n,l);
                    elseif (l==kappasz_pbk(s)+1)
                        M(snlj,end-ExNo) = M(snlj,end-ExNo) + bxy11snl(s,n,l); 
                        M(snlj,end-EyNo) = M(snlj,end-EyNo) + bxy21snl(s,n,l);
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bxy31snl(s,n,l);
                    end
                elseif (jj==l+1)
                        M(snlj,end-ExNo) = M(snlj,end-ExNo) + bxy12snl(s,n,l); 
                        M(snlj,end-EyNo) = M(snlj,end-EyNo) + bxy22snl(s,n,l);
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bxy32snl(s,n,l);
                end
            end
        end
    end
end

% Step 4: Jx or Jy
for s=1:S_pbk
    for ind_n=1:(2*Ns_pbk(s)+1)
        for l = 1:kappasz_pbk(s)+1
            snl1 = ind_M_pbk(s,ind_n,l,1);
            M(len_M_pbk,firstInd+snl1) = M(len_M_pbk,firstInd+snl1) + 1; 
        end
    end
end

% for x-dirction: ExyNo = ExNo
% M(len_subMxypbk,end-ExNo) = M(len_subMxypbk,end-ExNo) + bx10;
% for y-dirction:  ExyNo = EyNo
% M(len_subMxypbk,end-EyNo) = M(len_subMxypbk,end-EyNo) + by20;
M(len_M_pbk,end-ExyNo) = M(len_M_pbk,end-ExyNo) + bx10_by20;

end
