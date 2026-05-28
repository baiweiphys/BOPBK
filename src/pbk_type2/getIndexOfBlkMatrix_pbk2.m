function ind = getIndexOfBlkMatrix_pbk2(Ns_pbk,kappasz_pbk,MatrixNo)
% GETINDEXOFBLKMATRIX_PBK2 To obtain the index of a sub-matrix.
% @Description: To obtain the index of a sub-matrix within a composite 
% matrix that exhibits a PBK (type 2) plasma distribution.
% @Filename: getIndexOfBlkMatrix_pbk2.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-08-11
% @Modification: 
%   Bai Wei, 2026.04.04: update len_M_pbk2 for pbk type2.

% Calculate matrix size for PBK (type2) species.
% PBK type 2
len_M_pbk = 0.5*sum((2*Ns_pbk+1).*(kappasz_pbk+3).*kappasz_pbk) + 1;


switch MatrixNo
    case 1
        ind = 1:len_M_pbk;
    case 2
        ind = len_M_pbk+1:2*len_M_pbk;
    case 3
        ind = 2*len_M_pbk+1:3*len_M_pbk;   
    otherwise 
        disp("The PBK (type 2) MatrixNo must be an integer in the range 1 to 3.");
end


