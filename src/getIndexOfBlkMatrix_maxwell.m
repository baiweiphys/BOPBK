function idx = getIndexOfBlkMatrix_maxwell(S_bm,Ns_bm,J,MatrixNo)
% @Description: To obtain the index of a subblock matrix within a composite 
% matrix that exhibits a bi-Maxwellian plasma distribution.
% @Filename: getIndexOfBlkMatrix_maxwell.m
% @Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
% @Date: 2023-09-03
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-04

len_subMmaxwell = getLen_SNJmap2oneDim(S_bm,Ns_bm,J) + 1;

switch MatrixNo
    case 1
        idx = 1:len_subMmaxwell; % for Mx_maxwell
    case 2
        idx = len_subMmaxwell+1:2*len_subMmaxwell; % for My_maxwell
    case 3
        idx = 2*len_subMmaxwell+1:3*len_subMmaxwell; % for Mz_maxwell
    otherwise 
        disp('bi-Maxwellian MatrixNo must be an integer in the range 1 to 3');
end
