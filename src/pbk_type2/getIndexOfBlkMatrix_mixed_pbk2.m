function ind = getIndexOfBlkMatrix_mixed_pbk2(Ns_pbk,Ns_mk,kappasz_pbk,J,MatrixNo)
% GETINDEXOFBLKMATRIX_MIXED_PBK2 To obtain the index of a sub-matrix in the mixed solver. 
% @Description: To obtain the index of a subblock matrix within a composite 
% matrix that exhibits a mixed plasma distribution of product-bi-kappa (PBK type2) 
% and Maxwellian-kappa (MK).
% @Filename: getIndexOfBlkMatrix_mixed_pbk2.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-10-16
% @Modification: 
%   Bai Wei, 2026.04.02: Updated BM to MK.
%   Bai Wei, 2026.04.04: update len_M_pbk2 for pbk type2.


% PBK type 2
% Calculate matrix size for PBK (type2) species.
len_M_pbk = 0.5*sum((2*Ns_pbk+1).*(kappasz_pbk+3).*kappasz_pbk) + 1;

% Calculate matrix size for Maxwellian-Kappa (MK) species
len_M_mk = J*sum(2*Ns_mk+1) + 1;

switch MatrixNo
    case 1
        ind = 1:len_M_pbk; 
    case 2
        ind = len_M_pbk+1:2*len_M_pbk;
    case 3
        ind = 2*len_M_pbk+1:3*len_M_pbk;
    case 4
        ind = 3*len_M_pbk+1:3*len_M_pbk+len_M_mk;
    case 5
        ind = 3*len_M_pbk+len_M_mk+1:3*len_M_pbk+2*len_M_mk;
    case 6
        ind = 3*len_M_pbk+2*len_M_mk+1:3*len_M_pbk+3*len_M_mk;
    otherwise 
        disp("The mixed MatrixNo must be an integer in the range 1 to 6.");
end
