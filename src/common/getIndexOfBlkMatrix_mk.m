function ind = getIndexOfBlkMatrix_mk(Ns_mk,J,MatrixNo)
% GETINDEXOFBLKMATRIX_MK To obtain the index of a subblock matrix.
% @Description: To obtain the index of a subblock matrix within a composite 
% matrix that exhibits a Maxwellian-kappa (MK) plasma distribution.
% @Filename: getIndexOfBlkMatrix_mk.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-09-03
% @Modification: 
%   Bai Wei, 2026.04.02: Updated BM to MK.


% Calculate matrix size for bi-Maxwellian species
len_M_mk = J*sum(2*Ns_mk+1) + 1;

switch MatrixNo
    case 1
        ind = 1:len_M_mk; % for Mx_mk
    case 2
        ind = len_M_mk+1:2*len_M_mk; % for My_mk
    case 3
        ind = 2*len_M_mk+1:3*len_M_mk; % for Mz_mk
    otherwise 
        disp('MK MatrixNo must be a integrer for in range 1 to 3.');
end
