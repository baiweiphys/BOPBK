function M = M_mk(S_mk,Ns_mk,J,csnj,bxyz1snj,bxyz2snj,bxyz3snj, ...
                  MatrixNo,ExNo,EyNo,EzNo)
% M_MK Assemble the MK submatrix.
% @Description: Calculate matrix M_mk for the oblique 
% plasma wave model with a loss-cone bi-Maxwellian plasmas.
% @Filename: M_mk.m
% @Author: Bai Wei (baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.02: Updated to MK matrix. Structure identical to BM.

% common/PBK type 1 or 2
ind_M_mk = @(s,ind_n,j) SNJmap2oneDim_mk(Ns_mk,J,s,ind_n,j);

% Step 1
len_M_mk = J*sum(2*Ns_mk+1) + 1;

% Step 2: get the length of the Maxwellian-kappa (MK) matrix
% common/PBK type 1 or 2
ind_BlkMatrix = getIndexOfBlkMatrix_mk(Ns_mk,J,MatrixNo);
firstInd = ind_BlkMatrix(1)-1;

% Step 3: Create Matrix
len_row = len_M_mk;
len_col = 3*len_M_mk + 9;
M = zeros(len_row,len_col);

% Step 4: Assemble Matrix
for s=1:S_mk
    Nvec = -Ns_mk(s):Ns_mk(s);
    for ind_n=1:(2*Ns_mk(s)+1)
        n = Nvec(ind_n);
        for jj=1:J
            snj = ind_M_mk(s,ind_n,jj);
            % for xyz_{snj}
            M(snj,firstInd+snj) = M(snj,firstInd+snj) + csnj(s,n,jj);
            M(snj,end-ExNo) = M(snj,end-ExNo) + bxyz1snj(s,n,jj);
            M(snj,end-EyNo) = M(snj,end-EyNo) + bxyz2snj(s,n,jj);
            M(snj,end-EzNo) = M(snj,end-EzNo) + bxyz3snj(s,n,jj);
            % for Jxyz2
            M(len_M_mk,firstInd+snj) = M(len_M_mk,firstInd+snj) + csnj(s,n,jj);
            M(len_M_mk,end-ExNo) = M(len_M_mk,end-ExNo) + bxyz1snj(s,n,jj);
            M(len_M_mk,end-EyNo) = M(len_M_mk,end-EyNo) + bxyz2snj(s,n,jj);
            M(len_M_mk,end-EzNo) = M(len_M_mk,end-EzNo) + bxyz3snj(s,n,jj);
        end
    end
end
