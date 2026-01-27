function M = M_mixed_maxwell(S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,csnj, ...
                             bxyz1snj,bxyz2snj,bxyz3snj, ...
                             MatrixNo,ExNo,EyNo,EzNo)
% @Description: Calculate the M_mixed_maxwell matrix for the oblique plasma 
% wave model incorporating a mixed distribution of loss-cone PBK and 
% bi-Maxwellian plasmas.
% @Filename: M_mixed_maxwell.m
% @Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2026.01.25


idx_Mbm = @(s,idx_n,j) SNJmap2oneDim(Ns_bm,J,s,idx_n,j);

% Step 1
len_M_pbk = 0.5*sum((2*Ns_pbk+1).*(kappasz_pbk+4).*(kappasz_pbk+1)) + 1; % PBK submatrix
len_M_bm = J*sum(2*Ns_bm+1) + 1;

% Step 2: get the length of the maxwellian matrix
idx_BlkMatrix = getIndexOfBlkMatrix_mixed(Ns_pbk,Ns_bm,kappasz_pbk,J,MatrixNo);
firstIdx = idx_BlkMatrix(1)-1;


% Step 3: Create Matrix
% for PBK matrix 
len_pbk = 3*len_M_pbk;

% for bi-Maxwellian matrix 
len_bm = 3*len_M_bm;
%
len_row = len_M_bm;
len_col = len_pbk + len_bm + 9;
M = zeros(len_row,len_col);

% Step 4: Create Matrix
for s=1:S_bm
    Nvec = -Ns_bm(s):Ns_bm(s);
    for idx_n=1:(2*Ns_bm(s)+1)
        n = Nvec(idx_n);
        for jj=1:J
            snj = idx_Mbm(s,idx_n,jj);
            % for xyz_{snj}
            M(snj,firstIdx+snj) = M(snj,firstIdx+snj) + csnj(s,n,jj);
            M(snj,end-ExNo) = M(snj,end-ExNo) + bxyz1snj(s,n,jj); 
            M(snj,end-EyNo) = M(snj,end-EyNo) + bxyz2snj(s,n,jj);
            M(snj,end-EzNo) = M(snj,end-EzNo) + bxyz3snj(s,n,jj);
            % for Jxyz2
            M(len_M_bm,firstIdx+snj) = M(len_M_bm,firstIdx+snj) + csnj(s,n,jj);
            M(len_M_bm,end-ExNo) = M(len_M_bm,end-ExNo) + bxyz1snj(s,n,jj);
            M(len_M_bm,end-EyNo) = M(len_M_bm,end-EyNo) + bxyz2snj(s,n,jj);
            M(len_M_bm,end-EzNo) = M(len_M_bm,end-EzNo) + bxyz3snj(s,n,jj);
        end
    end
end
