function plasmaParams = getPlasmaParameters_pbk2(B0,par,J_opt) 
% GETPLASMAPARAMETERS_PBK2 To obtain the plasma parameters.
% @Description: To obtain the plasma parameters (PBK type 2).
% @Filename: getPlasmaParameters_pbk2.m
% @Date: 2022-01-05
% @Modification: 
%   Bai Wei, 2026.03.31: Added optional argument `is_pbk_vts` to Normalize 
%   using the PBK or bi-Maxwellian thermal velocity.
%   Bai Wei, 2026.04.02: Added optional argument `kappasx_th`.
%   Bai Wei, 2026.04.04: only update vtsz and vtsx for pbk type2.

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
    kappasz_threshold(s) = par.data(s,9);
    kappasx(s) = par.data(s,10);
    kappasx_threshold(s) = par.data(s,11);
    Ns(s) = par.data(s,12);
    is_pbk_vts(s) = par.data(s,13);
    As(s) = par.data(s,5)./par.data(s,4) - 1; % The s species temperature anisotropy
    
    if and(kappasz(s)<kappasz_threshold(s), is_pbk_vts(s))
        % if (kappa_sz(s) < kappa_threshold(s)) && is_pbk_vts(s)
        % thermal velocity of loss-cone PBK type 1 distribution
        % PBK type 2
        vtsz(s) = sqrt((2.0-3.0/kappasz(s))*kB*Tsz(s)/ms(s)); 
        vtsx(s) = sqrt((2.0-4.0/kappasx(s))*kB*Tsx(s)/ms(s)/(1.0+sgms(s)));
        % %
        % (M. S. dos Santos, 2024/2015/2026 PoP)
        % Normalize using the bi-Maxwellian thermal velocity
        % vtsz(s) = sqrt(2.0*kB*Tsz(s)/ms(s)); 
        % vtsx(s) = sqrt(2.0*kB*Tsx(s)/ms(s)/(1.0+sgms(s)));
    else
        % loss-cone bi-Maxwellian distribution
        vtsz(s) = sqrt(2.0*kB*Tsz(s)/ms(s)); 
        vtsx(s) = sqrt(2.0*kB*Tsx(s)/ms(s)/(1.0+sgms(s)));
    end
    % parallel drift velocity
    us0(s) = par.data(s,6).*sqrt(c2); 
end

% for i=1:S
%     % parallel drift velocity, normalized by vc
%     us0(s) = par.data(s,6).*sqrt(2*kB*Ts_parallel(3)/ms(3)); 
% end

lambdaDs = sqrt(epsilon0*kB*Tsz./(ns0.*qs.^2)); % Debye length, Tzs
kDs = 1.0./lambdaDs;

wps = sqrt(ns0.*qs.^2.0./ms/epsilon0); % plasma frequency
wcs = B0*qs./ms; % % cyclotron frequency
wce = abs(wcs(1));
rhocs = vtsx./abs(wcs); % cyclotron radius
wps2 = wps.^2;

is_pbk = (kappasz<kappasz_threshold); % index of PBK plasmas for s-th species
is_mk = (kappasz>kappasz_threshold); % index of Maxwellian-kappa (MK) plasmas for s-th species

% 
% betasz=2*mu0*kB.*ns0.*Tzs./B0^2;
% betasp=2*mu0*kB.*ns0.*Tps./B0^2;
% vA=B0/sqrt(mu0*sum(ms.*ns0)); %
% cS=sqrt(2*min(kB*Tzs)/max(ms));

%% Output parameters
plasmaParams.S = S;
plasmaParams.Ns = Ns;
plasmaParams.J_opt = J_opt;
plasmaParams.is_pbk = is_pbk;
plasmaParams.is_mk = is_mk;
plasmaParams.kappasz = kappasz;
plasmaParams.kappasx = kappasx;
plasmaParams.kappasx_threshold = kappasx_threshold; 
plasmaParams.vtsz = vtsz;
plasmaParams.vtsx = vtsx;
plasmaParams.Tsz = Tsz;
plasmaParams.Tsx = Tsx;
plasmaParams.sgms = sgms;
plasmaParams.wps = wps;
plasmaParams.wcs = wcs;
plasmaParams.us0 = us0;
plasmaParams.rhocs = rhocs;
plasmaParams.lambdaDs = lambdaDs;
plasmaParams.ms = ms;
plasmaParams.ns0 = ns0;
plasmaParams.is_pbk_vts = is_pbk_vts;

end
