function len = getSizeOfAllMatrix(Ns,J,index_pbk,index_bm,kappasz)
% @Description: To obtain the size of the complete matrix.
% @Filename: getSizeOfAllMatrix.m
% @Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
% @Date: 2024-09-20
% @LastEditors: Bai Wei
% @LastEditTime: 2025-02-01

S_pbk = sum(index_pbk);
S_bm = sum(index_bm);
Ns_pbk = Ns(index_pbk);
Ns_bm = Ns(index_bm);

kappasz_pbk = kappasz(index_pbk);

if S_bm==0
    tmp_pbk = getIndexOfBlkMatrix_pbk(S_pbk,Ns_pbk,kappasz_pbk,3);
    len = tmp_pbk(end) + 6;
elseif S_pbk==0
    tmp_bm = getIndexOfBlkMatrix_maxwell(S_bm,Ns_bm,J,3);
    len = tmp_bm(end) + 9;
else
    tmp_mixed = getIndexOfBlkMatrix_mixed(S_pbk,S_bm,Ns_pbk,Ns_bm,kappasz_pbk,J,6);
    len = tmp_mixed(end) + 9;
end






