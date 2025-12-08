function b5snl = pbk_b5snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,EPS0)
% @Description: Calculate the coefficients of b5snl for the oblique 
% plasma waves with a loss-cone PBK distribution.
% @Filename: pbk_b5snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-11-20
% @LastEditors: Bai Wei
% @LastEditTime: 2025-11-09

% FuncS_pbk = @(Jnum,dJnum,num,den) funcS_pbk(s,n,kappas_perp,lambdas,sgms,Jnum,dJnum,num,den,EPS0)
% S1_pbk = FuncS_pbk(2,0,1,2)
% S2_pbk = FuncS_pbk(2,0,1,1)
% S3_pbk = FuncS_pbk(1,1,2,2)
% S4_pbk = FuncS_pbk(1,1,2,1)
% S5_pbk = FuncS_pbk(0,2,3,2)
% S6_pbk = FuncS_pbk(0,2,3,1)
% S7_pbk = FuncS_pbk(2,0,-1,1)
% S8_pbk = FuncS_pbk(1,1,0,1)
% S9_pbk = FuncS_pbk(0,2,1,1)

S5 = @(s,n) funcS_pbk(s,n,kappasx,lambdas,sgms,0,2,3,2,EPS0);
S9 = @(s,n) funcS_pbk(s,n,kappasx,lambdas,sgms,0,2,1,1,EPS0);
S59 = @(s,n) (kappasx(s)+sgms(s)+1)/kappasx(s)*S5(s,n) ...
    -2*sgms(s)*lambdas(s)*S9(s,n);

bsnl = @(s,n,l) pbk_bsnl(s,n,l,kappasz,vtsz,wcs,kz);
b5snl = S59(s,n) * bsnl(s,n,l);

end