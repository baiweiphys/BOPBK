function b6snl = pbk1_b6snl(s,n,l,kz,kappasz,kappasx,vtsz,vtsx,sgms,lambdas,norm_omega,EPS0)
% PBK1_B6SNL Calculate the coefficient b6snl.
% @Description: Calculate the coefficients of b6snl for the oblique 
% plasma waves with a loss-cone PBK distribution.
% @Filename: pbk1_b6snl.m
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
S6 = @(s,n) funcS_pbk1(s,n,kappasx,lambdas,sgms,0,2,3,1,EPS0);

% PBK type 1
% bsl = @(s,l) pbk1_bsl(s,l,kappasz,vtsz,vtsx,kz,norm_omega);
b6snl = S6(s,n) * pbk1_bsl(s,l,kappasz,vtsz,vtsx,kz,norm_omega);

end
