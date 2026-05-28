function M = M_mk_mixed_pbk2(S_mk,Ns_pbk,Ns_mk,kappasz_pbk,J,csnj, ...
                        bxyz1snj,bxyz2snj,bxyz3snj, ...
                        MatrixNo,ExNo,EyNo,EzNo)
% M_MK_MIXED_PBK2 Assemble the MK submatrix in the mixed solver.
% @Description: Calculate the M_mixed_mk matrix for the oblique plasma 
% wave model incorporating a mixed distribution of loss-cone PBK (type 2) 
% and Maxwellian-kappa (MK) plasmas.
% @Filename: M_mk_mixed_pbk2.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.02: Updated to MK matrix. Structure identical to BM.
%   Bai Wei, 2026.04.04: update len_M_pbk2 for pbk type2.


% common/PBK type 1 or 2
ind_M_mk = @(s,ind_n,j) SNJmap2oneDim_mk(Ns_mk,J,s,ind_n,j);

% Step 1
% PBK type 2
len_M_pbk = 0.5*sum((2*Ns_pbk+1).*(kappasz_pbk+3).*kappasz_pbk) + 1; % PBK2 submatrix
len_M_mk = J*sum(2*Ns_mk+1) + 1;

% Step 2: get the length of the mk matrix
% PBK type 2
ind_BlkMatrix = getIndexOfBlkMatrix_mixed_pbk2(Ns_pbk,Ns_mk,kappasz_pbk,J,MatrixNo);
firstInd = ind_BlkMatrix(1)-1;


% Step 3: Create Matrix
% for PBK2 matrix 
len_pbk = 3*len_M_pbk;

% for Maxwellian-kappa (MK) matrix 
len_mk = 3*len_M_mk;
%
len_row = len_M_mk;
len_col = len_pbk + len_mk + 9;
M = zeros(len_row,len_col);

% Step 4: Create Matrix
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
