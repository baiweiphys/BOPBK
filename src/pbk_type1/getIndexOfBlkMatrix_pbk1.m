function ind = getIndexOfBlkMatrix_pbk1(Ns_pbk,kappasz_pbk,MatrixNo)
% GETINDEXOFBLKMATRIX_PBK1 To obtain the index of a sub-matrix.
% @Description: To obtain the index of a sub-matrix within a composite 
% matrix that exhibits a PBK (type 1) plasma distribution.
% @Filename: getIndexOfBlkMatrix_pbk1.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-08-11
% @Modification: 
%   Bai Wei, 2026.04.05: Updated to PBK type 1.


% Calculate matrix size for PBK (type 1) species
% PBK type 1
len_M_pbk = 0.5*sum((2*Ns_pbk+1).*(kappasz_pbk+4).*(kappasz_pbk+1)) + 1;

switch MatrixNo
    case 1
        ind = 1:len_M_pbk;
    case 2
        ind = len_M_pbk+1:2*len_M_pbk;
    case 3
        ind = 2*len_M_pbk+1:3*len_M_pbk;   
    otherwise 
        disp("The PBK (type 1) MatrixNo must be an integer in the range 1 to 3.");
end


