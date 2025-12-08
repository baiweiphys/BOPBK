function M = Mz_pbk(S_pbk,Ns_pbk,kappasz_pbk,csn_pbk, ...
                    bz11snl,bz21snl,bz31snl, ...
                    bz13snl,bz23snl,bz33snl, ...
                    bz12snl,bz22snl,bz32snl, ...
                    bz34snl,by20,MatrixNo,ExNo,EyNo,EzNo)
% @Description: To obtain matrix Mz_pbk for loss-cone PBK plasmas 
% in z-direction.
% @Filename: Mz_pbk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-08-12
% @LastEditors: Bai Wei
% @LastEditTime: 2025.02.13

index_Mpbk = @(s,idx_n,l,j) SNLJmap2oneDim(Ns_pbk,kappasz_pbk,s,idx_n,l,j);
len_Mpbk = getLen_SNLJmap2oneDim(S_pbk,Ns_pbk,kappasz_pbk) + 1;
index_BlkMatrix = getIndexOfBlkMatrix_pbk(S_pbk,Ns_pbk,kappasz_pbk,MatrixNo);
firstIndex = index_BlkMatrix(1)-1;
len_row = len_Mpbk;
len_col = 3*len_Mpbk + 6;
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
                    M(snlj,firstIndex+snlj+1) = M(snlj,firstIndex+snlj+1)+1;
                end
                if (jj==l)
                    if (l<=kappasz_pbk(s)-1)
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bz34snl(s,n,l+2); 
                        M(snlj,end-ExNo) = M(snlj,end-ExNo) + bz13snl(s,n,l+1); 
                        M(snlj,end-EyNo) = M(snlj,end-EyNo) + bz23snl(s,n,l+1);
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bz33snl(s,n,l+1);
                        M(snlj,end-ExNo) = M(snlj,end-ExNo) + bz11snl(s,n,l);
                        M(snlj,end-EyNo) = M(snlj,end-EyNo) + bz21snl(s,n,l);
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bz31snl(s,n,l);
                    elseif (l==kappasz_pbk(s))
                        M(snlj,end-ExNo) = M(snlj,end-ExNo) + bz13snl(s,n,l+1); 
                        M(snlj,end-EyNo) = M(snlj,end-EyNo) + bz23snl(s,n,l+1); 
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bz33snl(s,n,l+1);
                        M(snlj,end-ExNo) = M(snlj,end-ExNo) + bz11snl(s,n,l);
                        M(snlj,end-EyNo) = M(snlj,end-EyNo) + bz21snl(s,n,l);
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bz31snl(s,n,l);
                    elseif (l==kappasz_pbk(s)+1)
                        M(snlj,end-ExNo) = M(snlj,end-ExNo) + bz11snl(s,n,l);
                        M(snlj,end-EyNo) = M(snlj,end-EyNo) + bz21snl(s,n,l);
                        M(snlj,end-EzNo) = M(snlj,end-EzNo) + bz31snl(s,n,l);
                    end
                elseif (jj==l+1)
                    M(snlj,end-ExNo) = M(snlj,end-ExNo) + bz12snl(s,n,l); 
                    M(snlj,end-EyNo) = M(snlj,end-EyNo) + bz22snl(s,n,l); 
                    M(snlj,end-EzNo) = M(snlj,end-EzNo) + bz32snl(s,n,l); 
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

M(len_Mpbk,end-EzNo) = M(len_Mpbk,end-EzNo) + by20;

end
