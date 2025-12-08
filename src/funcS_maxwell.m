function funcS = funcS_maxwell(s,n,lambdas,sgms,Jnum,dJnum,num,EPS0)
% @Description: Compute the integral of funcS_maxwell for maxwellian of loss-cone PBK model.
% @Filename: funcS_maxwell.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-11-20
% @LastEditors: Bai Wei
% @LastEditTime: 2025-01-20

% example:
% funcS = @(Jnum,dJnum,num) funcS_maxwell(s,n,lambdas,sgms,Jnum,dJnum,num,EPS0)
% S1_maxwell = S2_maxwell = funcS(2,0,1)
% S3_maxwell = S4_maxwell = funcS(1,1,2)
% S5_maxwell = S6_maxwell = funcS(0,2,3)
% S7_maxwell = funcS(2,0,-1)
% S8_maxwell = funcS(1,1,0)
% S9_maxwell = funcS(0,2,1)

lambda = lambdas(s);
sgm = sgms(s);

Jn = @(x) besselj(n,x);
dJn = @(x) -0.5*besselj(n+1,x) + 0.5*besselj(n-1,x);

S0 = 4.0/(2*lambda)^(sgm+2)/gamma(sgm+1);

xmin = 0;
% xmax = 40;
xmax = Inf;
%
reltol=1e-16; 
abstol=1e-14;

eps0 = 1e-10;
% eps0 = eps;

if lambda > EPS0
    tmp = @(x) (x+eps0).^(2*sgm+num).*Jn(x).^Jnum.*dJn(x).^dJnum.*exp(-0.5*x.^2/lambda);
    % tmp_int = integral(tmp,xmin,xmax);
    tmp_int = integral(tmp,xmin,xmax,'RelTol',reltol,'AbsTol',abstol);
else
    tmp = @(x) (x*sqrt(2*lambda)+eps0).^(2*sgm+num)*sqrt(2*lambda)....
        .*Jn(x*sqrt(2*lambda)).^Jnum.*dJn(x*sqrt(2*lambda)).^dJnum.*exp(-x.^2);
    % tmp_int = integral(tmp,xmin,xmax/sqrt(2*lambda));
    tmp_int = integral(tmp,xmin,xmax/sqrt(2*lambda),'RelTol',reltol,'AbsTol',abstol);
end

funcS = S0*tmp_int;

end