function [dy_logdt,Gamma,Omega,rr,Sph,kr] = DaeFunc(t,y_log, ...
    Kelet,ExtraRho,K_e,Add_e,np,ns,v,...
    fAssembleKD,fAssembleKMu,fApplyBC,fComputeTe,fComputeCoefficients,fComputeOmega,fComputePhotoSource,...
    Gamma2GammaClose,N2RhoEps,BCel_val,indices_boundaries)
%MYODEFUNC Summary of this function goes here
%   Detailed explanation goes here

% y = exp(y_log);
y = y_log;

n = y(1:np*ns,:);
phi = y(np*ns+1:end,:);

% The block of code below should be the body of the function "MyOdeFunc"
% except for the computation of E that in the DAE case should be
% E = K_e * phi + Add_e * BCel_val(t);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% START
rho_eps = N2RhoEps * n;

E = K_e * phi + Add_e * (BCel_val*v(t));

Te = fComputeTe(E);

[kr, mu, D] = fComputeCoefficients(E,Te);

KD = fAssembleKD(D(:));
[KMu,u] = fAssembleKMu(E,mu(:));

Gamma = (KD + KMu) * n;
[~,rr,s] = fComputeOmega(n,kr,E);
Sph = fComputePhotoSource(E,s);

Omega = reshape(fComputeOmega(n, kr, E),[],1);

rr = reshape(rr,[],1);
kr = reshape(kr,[],1);

% apply BC
GammaBC = fApplyBC(n(indices_boundaries.nodes), ...
                   u(indices_boundaries.interfaces), ...
                   Te,t,E);
Gamma(indices_boundaries.interfaces) = GammaBC;

GammaClose = Gamma2GammaClose * Gamma;

dndt = -GammaClose + Omega;
% END
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

dy_logdt = [dndt; Kelet*phi - rho_eps - ExtraRho*(BCel_val*v(t))];

end
