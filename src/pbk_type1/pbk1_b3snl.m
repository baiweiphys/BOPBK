function b3snl = pbk1_b3snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,norm_omega,EPS0)
% PBK1_B3SNL Calculate the coefficient b3snl. 
% @Description: Calculate the coefficients of b3snl for the oblique 
% plasma waves with a loss-cone PBK distribution.
% @Filename: pbk1_b3snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-11-20
% @Modification: 
%   Bai Wei, 2026.04.01: add norm_omega.


%% FuncS_pbk = @(Jnum,dJnum,num,den) funcS_pbk1(s,n,kappas_perp,lambdas,sgms,Jnum,dJnum,num,den,pbk_type,EPS0)
% S1_pbk = FuncS_pbk1(2,0,1,2)
% S2_pbk = FuncS_pbk1(2,0,1,1)
% S3_pbk = FuncS_pbk1(1,1,2,2)
% S4_pbk = FuncS_pbk1(1,1,2,1)
% S5_pbk = FuncS_pbk1(0,2,3,2)
% S6_pbk = FuncS_pbk1(0,2,3,1)
% S7_pbk = FuncS_pbk1(2,0,-1,1)
% S8_pbk = FuncS_pbk1(1,1,0,1)
% S9_pbk = FuncS_pbk1(0,2,1,1)

% PBK type 1
S3 = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,1,1,2,2,EPS0);
S8 = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,1,1,0,1,EPS0);
S38 = @(s,n) (kappasx(s)+sgms(s)+1)/kappasx(s)*S3(s,n)...
    -2*sgms(s)*lambdas(s)*S8(s,n);

% PBK type 1
% bsnl = @(s,n,l) pbk1_bsnl(s,n,l,kappasz,vtsz,wcs,kz,norm_omega);
b3snl = S38(s,n) * pbk1_bsnl(s,n,l,kappasz,vtsz,wcs,kz,norm_omega);

end
