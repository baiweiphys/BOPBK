function M = M_mixed_maxwell(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csnj, ...
                             bxyz1snj,bxyz2snj,bxyz3snj, ...
                             MatrixNo,ExNo,EyNo,EzNo)
% @Description: Calculate the M_mixed_maxwell matrix for the oblique plasma 
% wave model incorporating a mixed distribution of loss-cone PBK and 
% bi-Maxwellian plasmas.
% @Filename: M_mixed_maxwell.m
% @Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025.02.13

index_Mbm = @(s,idx_n,j) SNJmap2oneDim(Ns_bm,J,s,idx_n,j);

len_Mpbk = getLen_SNLJmap2oneDim(S_pbk,Ns_pbk,kappasz_pbk)+1;
len_Mbm = getLen_SNJmap2oneDim(S_bm,Ns_bm,J) + 1;
index_BlkMatrix = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,MatrixNo);
firstIndex = index_BlkMatrix(1)-1;
len_pbk = 3*len_Mpbk;
len_bm = 3*len_Mbm;
len_row = len_Mbm;
len_col = len_pbk + len_bm + 9;
M = zeros(len_row,len_col);
for s=1:S_bm
    Nvec = -Ns_bm(s):Ns_bm(s);
    for idx_n=1:(2*Ns_bm(s)+1)
        n = Nvec(idx_n);
        for jj=1:J
            snj = index_Mbm(s,idx_n,jj);
            M(snj,firstIndex+snj) = M(snj,firstIndex+snj) + csnj(s,n,jj);
            M(snj,end-ExNo) = M(snj,end-ExNo) + bxyz1snj(s,n,jj); 
            M(snj,end-EyNo) = M(snj,end-EyNo) + bxyz2snj(s,n,jj);
            M(snj,end-EzNo) = M(snj,end-EzNo) + bxyz3snj(s,n,jj);
            M(len_Mbm,firstIndex+snj) = M(len_Mbm,firstIndex+snj) + csnj(s,n,jj);
            M(len_Mbm,end-ExNo) = M(len_Mbm,end-ExNo) + bxyz1snj(s,n,jj);
            M(len_Mbm,end-EyNo) = M(len_Mbm,end-EyNo) + bxyz2snj(s,n,jj);
            M(len_Mbm,end-EzNo) = M(len_Mbm,end-EzNo) + bxyz3snj(s,n,jj);
        end
    end
end
