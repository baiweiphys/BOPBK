function funcS = funcS_pbk2(s,n,kappasx,lambdas,sgms,p1,p2,p3,p4,EPS0)
% FUNCS_PBK2 Compute the integral of funcS_pbk2.
% @Description: Compute the integral of funcS_pbk for loss-cone PBK model (type 2).
% @Filename: funcS_pbk2.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-11-20
% @Modification: 
%   Bai Wei, 2026.03.31: FuncS_pbk Optimization.
%   Bai Wei, 2026.04.04: only update S0 for pbk type2.


%% example:
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

%%
lambda = lambdas(s);
kappax = kappasx(s);
sgm = sgms(s);

Jn = @(x) besselj(n,x);
dJn = @(x) -0.5*besselj(n+1,x) + 0.5*besselj(n-1,x);

xmin = 0;
% xmax = 30;
xmax = Inf;
%
% reltol=1e-14; 
% abstol=1e-10;

eps0 = eps; 
% eps0 = 1.0e-14;


if lambda > EPS0
    % PBK type 2
    S0 = 4.0*(2*lambda)^(-sgm-2)*kappax^(-sgm-1) ...
        *exp(gammaln(kappax+sgm)-gammaln(kappax-1)) / gamma(sgm+1);
    %
    integrand = @(x) (x+eps0).^(2*sgm+p3).*Jn(x).^p1.*dJn(x).^p2...
            ./((1+0.5*x.^2/lambda/kappax).^(kappax+sgm+p4));
    tmp_int = integral(integrand,xmin,xmax);
    % tmp_int = integral(tmp,xmin,xmax,'RelTol',reltol,'AbsTol',abstol);
    % tmp_int = quadgk(tmp,xmin,xmax);
    % tmp_int = quadgk(tmp,xmin,xmax,'MaxIntervalCount', 2000);
    funcS = S0*tmp_int;
else
    % PBK type 2
    S0 = 4.0*(2*lambda)^(p3/2.0-1.5)*kappax^(-sgm-1) ...
        *exp(gammaln(kappax+sgm)-gammaln(kappax-1)) / gamma(sgm+1);
    %
    integrand = @(x) (x+eps0).^(2*sgm+p3) ...
            .*Jn(x*sqrt(2*lambda)).^p1.*dJn(x*sqrt(2*lambda)).^p2 ...
            ./((1+x.^2/kappax).^(kappax+sgm+p4));
    tmp_int = integral(integrand,xmin,xmax/sqrt(2*lambda));
    % tmp_int = integral(tmp,xmin,xmax/sqrt(2*lambda),'RelTol',reltol,'AbsTol',abstol);
    % tmp_int = quadgk(tmp,xmin,xmax/sqrt(2*lambda));
    % tmp_int = quadgk(tmp,xmin,xmax/sqrt(2*lambda),'MaxIntervalCount', 2000);
    funcS = S0*tmp_int;
end
