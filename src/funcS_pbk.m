function funcS = funcS_pbk(s,n,kappasx,lambdas,sgms,p1,p2,p3,p4,EPS0)
% @Description: Compute the integral of funcS_pbk for loss-cone PBK model.
% @Filename: funcS_pbk.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-11-20
% @LastEditors: Bai Wei
% @LastEditTime: 2025-11-09

% example:
% funcS = @(Jnum,dJnum,num,den) funcS_pbk(s,n,kappas_perp,lambdas,sgms,Jnum,dJnum,num,den,EPS0)
% S1_pbk = FuncS_pbk(2,0,1,2)
% S2_pbk = FuncS_pbk(2,0,1,1)
% S3_pbk = FuncS_pbk(1,1,2,2)
% S4_pbk = FuncS_pbk(1,1,2,1)
% S5_pbk = FuncS_pbk(0,2,3,2)
% S6_pbk = FuncS_pbk(0,2,3,1)
% S7_pbk = FuncS_pbk(2,0,-1,1)
% S8_pbk = FuncS_pbk(1,1,0,1)
% S9_pbk = FuncS_pbk(0,2,1,1)

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
    S0 = 4.0*(2*lambda)^(-sgm-2)*kappax^(-sgm-1)...
    *exp(gammaln(kappax+sgm+1)-gammaln(kappax)) / gamma(sgm+1);
    %
    integrand = @(x) (x+eps0).^(2*sgm+p3).*Jn(x).^p1.*dJn(x).^p2...
            ./((1+0.5*x.^2/lambda/kappax).^(kappax+sgm+p4));
    tmp_int = integral(integrand,xmin,xmax);
    % tmp_int = integral(tmp,xmin,xmax,'RelTol',reltol,'AbsTol',abstol);
    % tmp_int = quadgk(tmp,xmin,xmax);
    % tmp_int = quadgk(tmp,xmin,xmax,'MaxIntervalCount', 2000);
    funcS = S0*tmp_int;
else
    S0 = 4.0*(2*lambda)^(p3/2.0-1.5)*kappax^(-sgm-1) ...
    *exp(gammaln(kappax+sgm+1)-gammaln(kappax)) / gamma(sgm+1);
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
