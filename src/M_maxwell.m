function M = M_maxwell(S_bm,Ns_bm,J,csnj,bxyz1snj,bxyz2snj,bxyz3snj, ...
                       MatrixNo,ExNo,EyNo,EzNo)
% @Description: Calculate matrix M_maxwell for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian plasmas.
% @Filename: M_maxwell.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-13


idx_M_bm = @(s,idx_n,j) SNJmap2oneDim(Ns_bm,J,s,idx_n,j);


% Step 1
len_M_bm = J*sum(2*Ns_bm+1) + 1;

% Step 2: get the length of the maxwellian matrix
index_BlkMatrix = getIndexOfBlkMatrix_maxwell(Ns_bm,J,MatrixNo);
firstIndex = index_BlkMatrix(1)-1;


% Step 3: Create Matrix
len_row = len_M_bm;
len_col = 3*len_M_bm + 9;
M = zeros(len_row,len_col);


% Step 4: Assemble Matrix
for s=1:S_bm
    Nvec = -Ns_bm(s):Ns_bm(s);
    for idx_n=1:(2*Ns_bm(s)+1)
        n = Nvec(idx_n);
        for jj=1:J
            snj = idx_M_bm(s,idx_n,jj);
            % for xyz_{snj}
            M(snj,firstIndex+snj) = M(snj,firstIndex+snj) + csnj(s,n,jj);
            M(snj,end-ExNo) = M(snj,end-ExNo) + bxyz1snj(s,n,jj);
            M(snj,end-EyNo) = M(snj,end-EyNo) + bxyz2snj(s,n,jj);
            M(snj,end-EzNo) = M(snj,end-EzNo) + bxyz3snj(s,n,jj);
            % for Jxyz2
            M(len_M_bm,firstIndex+snj) = M(len_M_bm,firstIndex+snj) + csnj(s,n,jj);
            M(len_M_bm,end-ExNo) = M(len_M_bm,end-ExNo) + bxyz1snj(s,n,jj);
            M(len_M_bm,end-EyNo) = M(len_M_bm,end-EyNo) + bxyz2snj(s,n,jj);
            M(len_M_bm,end-EzNo) = M(len_M_bm,end-EzNo) + bxyz3snj(s,n,jj);
        end
    end
end
