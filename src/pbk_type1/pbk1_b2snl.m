function b2snl = pbk1_b2snl(s,n,l,kz,kappasz,kappasx,vtsz,vtsx,sgms,lambdas,norm_omega,EPS0)
% PBK1_B2SNL Calculate the coefficient b2snl.
% @Description: Calculate the coefficients of b2snl for the oblique 
% plasma waves with a loss-cone PBK distribution.
% @Filename: pbk1_b2snl.m
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
S2 = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,2,0,1,1,EPS0);

% PBK type 1
% bsl = @(s,l) pbk1_bsl(s,l,kappasz,vtsz,vtsx,kz,norm_omega);
b2snl = S2(s,n) * pbk1_bsl(s,l,kappasz,vtsz,vtsx,kz,norm_omega);

end
