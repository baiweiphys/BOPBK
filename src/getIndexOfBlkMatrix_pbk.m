function idx = getIndexOfBlkMatrix_pbk(S_pbk,Ns_pbk,kappasz_pbk,MatrixNo)
% @Description: To obtain the index of a subblock matrix within a composite 
% matrix that exhibits a PBK plasma distribution.
% @Filename: getIndexOfBlkMatrix_pbk.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-08-11
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-01

len_subMpbk = getLen_SNLJmap2oneDim(S_pbk,Ns_pbk,kappasz_pbk) + 1;

switch MatrixNo
    case 1
        idx = 1:len_subMpbk;
    case 2
        idx = len_subMpbk+1:2*len_subMpbk;
    case 3
        idx = 2*len_subMpbk+1:3*len_subMpbk;   
    otherwise 
        disp("The PBK MatrixNo must be an integer in the range 1 to 3.");
end


