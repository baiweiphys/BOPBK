function [S,Ns,index_pbk,index_bm,kappasz,kappasx,vtsz,vtsx,Tsz,Tsx,sgms,wps,wcs,us0,rhocs,lambdaDs] = getPlasmaParameters(B0,par) 
% @Description: To obtain the plasma parameters.
% @Filename: getPlasmaParameters.m
% @Date: 2022-01-05
% @LastEditors: Bai Wei
% @LastEditTime: 2025-11-07

params_with_unit;

[S,~] = size(par.data);

for s=1:S
    qs(s)=par.data(s,1)*qe; % charge
    ms(s)=par.data(s,2)*me; % mass
    ns0(s)=par.data(s,3); % desity unit: m^-3
    Tsz(s) = par.data(s,4)*qe/kB; % parallel temperature, unit: eV -> K
    Tsx(s) = par.data(s,5)*qe/kB; % perp temperature, unit: eV -> K
    sgms(s) = par.data(s,7); % loss-cone: sigma
    kappasz(s) = par.data(s,8);
    kappasx(s) = par.data(s,9);
    kappasz_threshold(s) = par.data(s,10);
    Ns(s) = par.data(s,11);
    As(s) = par.data(s,5)./par.data(s,4) - 1; % The s species temperature anisotropy
    
    if (kappasz(s)<kappasz_threshold(s))
        % thermal velocity of loss-cone PBK distribution
        vtsz(s) = sqrt((2.0-1.0/kappasz(s))*kB*Tsz(s)/ms(s)); 
        vtsx(s) = sqrt((2.0-2.0/kappasx(s))*kB*Tsx(s)/ms(s)/(1.0+sgms(s)));
    else
        % loss-cone bi-Maxwellian distribution
        vtsz(s) = sqrt(2.0*kB*Tsz(s)/ms(s)); 
        vtsx(s) = sqrt(2.0*kB*Tsx(s)/ms(s)/(1.0+sgms(s)));
    end
    %
    % (M. S. dos Santos, 2016 PoP)
    % for loss-cone PBK and bi-Maxwellian distribution
    % vtsz(s) = sqrt(2.0*kB*Tsz(s)/ms(s)); 
    % vtsx(s) = sqrt(2.0*kB*Tsx(s)/ms(s)/(1.0+sgms(s)));
    %
    % parallel drift velocity
    us0(s) = par.data(s,6).*sqrt(c2); 
end


lambdaDs = sqrt(epsilon0*kB*Tsz./(ns0.*qs.^2)); % Debye length, Tzs
kDs = 1.0./lambdaDs;

wps = sqrt(ns0.*qs.^2.0./ms/epsilon0); % plasma frequency
wcs = B0*qs./ms; % % cyclotron frequency
wce = abs(wcs(1));
rhocs = vtsx./abs(wcs); % cyclotron radius
wps2 = wps.^2;

index_pbk = (kappasz<kappasz_threshold); % index of PBK plasmas for s-th species
index_bm = (kappasz>kappasz_threshold); % index of bi-Maxwellian plasmas for s-th species


end