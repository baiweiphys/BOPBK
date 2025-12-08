function M = Mxy_pbk_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csn_pbk, ...
                           bxy11snl,bxy21snl,bxy31snl, ...
                           bxy12snl,bxy22snl,bxy32snl, ...
                           bxy33snl,bx10_by20, ...
                           MatrixNo,ExNo,EyNo,EzNo,ExyNo)
% @Description: To obtain matrix Mxy_pbk_mixed for the mixed loss-cone
% distributions of PBK and bi-Maxwellian plasmas in x- or y-direction.
% @Filename: Mxy_pbk_mixed.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-10-16
% @LastEditors: Bai Wei
% @LastEditTime: 2025.02.13

index_Mpbk = @(s,idx_n,l,j) SNLJmap2oneDim(Ns_pbk,kappasz_pbk,s,idx_n,l,j);
len_Mpbk = getLen_SNLJmap2oneDim(S_pbk,Ns_pbk,kappasz_pbk)+1;
len_Mbm = getLen_SNJmap2oneDim(S_bm,Ns_bm,J)+1; 

index_BlkMatrix = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,MatrixNo);
firstIndex = index_BlkMatrix(1)-1;
len_pbk = 3*len_Mpbk;
len_bm = 3*len_Mbm;
len_row = len_Mpbk;
len_col = len_pbk + len_bm + 9;
M = zeros(len_row,len_col);
for s=1:S_pbk
    Nvec = -Ns_pbk(s):Ns_pbk(s);
    for idx_n=1:(2*Ns_pbk(s)+1)
        n = Nvec(idx_n);
        for l = 1:kappasz_pbk(s)+1
            for jj = 1:l+1
                snlj = index_Mpbk(s,idx_n,l,jj);
                M(snlj,firstIndex+snlj) = M(snlj,firstIndex+snlj) + csn_pbk(s,n); 
                if(jj<l+1) 
                    %snljp1 = index_Mpbk(s,idx_n,l,jj+1);
                    M(snlj,firstIndex+snlj+1) = M(snlj,firstIndex+snlj+1)+1;
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


for s=1:S_pbk
    for idx_n=1:(2*Ns_pbk(s)+1)
        for l = 1:kappasz_pbk(s)+1
            snl1 = index_Mpbk(s,idx_n,l,1);
            M(len_Mpbk,firstIndex+snl1) = M(len_Mpbk,firstIndex+snl1) + 1; 
        end
    end
end
M(len_Mpbk,end-ExyNo) = M(len_Mpbk,end-ExyNo) + bx10_by20;

end
