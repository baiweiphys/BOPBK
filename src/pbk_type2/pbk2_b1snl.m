function b1snl = pbk2_b1snl(s,n,l,kz,kappasz,kappasx,vtsz,sgms,wcs,lambdas,norm_omega,EPS0)
% PBK2_B1SNL Calculate the coefficient b1snl.
% @Description: Calculate the coefficients of b1snl for the oblique 
% plasma waves with a loss-cone PBK distribution (type 2).
% @Filename: pbk2_b1snl.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-11-20
% @Modification: 
%   Bai Wei, 2026.04.01: add norm_omega.
%   Bai Wei, 2026.04.04: only update S1, S2, and S17 for pbk type2.

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

% PBK type 2
S1 = @(s,n) funcS_pbk2(s,n,kappasx,lambdas,sgms,2,0,1,1,EPS0);
S7 = @(s,n) funcS_pbk2(s,n,kappasx,lambdas,sgms,2,0,-1,0,EPS0);
S17 = @(s,n) (kappasx(s)+sgms(s))/kappasx(s)*S1(s,n) ...
        -2*sgms(s)*lambdas(s)*S7(s,n);

% PBK type 2
% bsnl = @(s,n,l) pbk2_bsnl(s,n,l,kappasz,vtsz,wcs,kz,norm_omega);
b1snl = S17(s,n) * pbk2_bsnl(s,n,l,kappasz,vtsz,wcs,kz,norm_omega);

end
