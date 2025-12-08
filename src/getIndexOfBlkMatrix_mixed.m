function idx = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,MatrixNo)
% @Description: To obtain the index of a subblock matrix within a composite 
% matrix that exhibits a mixed plasma distribution of PBK and bi-Maxwellian.
% @Filename: getIndexOfBlkMatrix_mixed.m
% @Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
% @Date: 2023-10-16
% @LastEditors: Bai Wei
% @LastEditTime: 2025.02.13

len_Mpbk = getLen_SNLJmap2oneDim(S_pbk,Ns_pbk,kappasz_pbk) + 1;
len_Mmaxwell = getLen_SNJmap2oneDim(S_bm,Ns_bm,J) + 1;

switch MatrixNo
    case 1
        idx = 1:len_Mpbk; 
    case 2
        idx = len_Mpbk+1:2*len_Mpbk;
    case 3
        idx = 2*len_Mpbk+1:3*len_Mpbk;
    case 4
        idx = 3*len_Mpbk+1:3*len_Mpbk+len_Mmaxwell;
    case 5
        idx = 3*len_Mpbk+len_Mmaxwell+1:3*len_Mpbk+2*len_Mmaxwell;
    case 6
        idx = 3*len_Mpbk+2*len_Mmaxwell+1:3*len_Mpbk+3*len_Mmaxwell;
    otherwise 
        disp("The mixed MatrixNo must be an integer in the range 1 to 6.");
end
