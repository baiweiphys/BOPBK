function b5snl = pbk2_b5snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,norm_omega,EPS0)
% PBK2_B5SNL Calculate the coefficient b5snl.
% @Description: Calculate the coefficients of b5snl for the oblique 
% plasma waves with a PBK distribution (type 2).
% @Filename: pbk2_b5snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-11-20
% @Modification: 
%   Bai Wei, 2026.04.01: add norm_omega.
%   Bai Wei, 2026.04.04: only update S5, S9, and S59 for pbk type2.


% FuncS_pbk = @(Jnum,dJnum,num,den) funcS_pbk2(s,n,kappas_perp,lambdas,sgms,Jnum,dJnum,num,den,EPS0)
% S1_pbk = FuncS_pbk2(2,0,1,1)
% S2_pbk = FuncS_pbk2(2,0,1,0)
% S3_pbk = FuncS_pbk2(1,1,2,1)
% S4_pbk = FuncS_pbk2(1,1,2,0)
% S5_pbk = FuncS_pbk2(0,2,3,1)
% S6_pbk = FuncS_pbk2(0,2,3,0)
% S7_pbk = FuncS_pbk2(2,0,-1,0)
% S8_pbk = FuncS_pbk2(1,1,0,0)
% S9_pbk = FuncS_pbk2(0,2,1,0)

% PBK type 2% PBK type 2
S5 = @(s,n) funcS_pbk2(s,n,kappasx,lambdas,sgms,0,2,3,1,EPS0);
S9 = @(s,n) funcS_pbk2(s,n,kappasx,lambdas,sgms,0,2,1,0,EPS0);
S59 = @(s,n) (kappasx(s)+sgms(s))/kappasx(s)*S5(s,n) ...
    -2*sgms(s)*lambdas(s)*S9(s,n);


% PBK type 2
% bsnl = @(s,n,l) pbk2_bsnl(s,n,l,kappasz,vtsz,wcs,kz,norm_omega);
b5snl = S59(s,n) * pbk2_bsnl(s,n,l,kappasz,vtsz,wcs,kz,norm_omega);

end
