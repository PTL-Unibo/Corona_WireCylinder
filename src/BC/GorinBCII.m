function [GammaBC] = GorinBCII(n_boundary,u_boundary,r,gamma,qs,v_th_single_line,ns,Te)
%APPLYGORINBC Summary of this function goes here
%   Detailed explanation goes here
kB = 1.380649e-23;
me = 9.1093837e-31;
% Taken from "Boundary conditions for drift-diffusion equations in gas-discharge plasmas"
% Equations (9) and (10)
% DOI: https://doi.org/10.1063/1.5120613
n_boundary = reshape(n_boundary,2,ns);
u_boundary = reshape(u_boundary,2,ns);

u_boundary(1,:) = -u_boundary(1,:);
% the velocity u_boundary is pointing towards the wall

Te = Te*11600;
v_th_e = sqrt(8*kB*Te(1)./(pi*me));
v_th_single_line(1,1) = v_th_e;

v_th = repmat(v_th_single_line,2,1);
v_E = max(u_boundary,0);

GammaBC = (1-r)./(1+r) .* n_boundary .* (0.5*v_th + v_E);

Gamma_plus_collettore = gamma * sum(GammaBC(end,qs>0),2);
GammaBC(end,1) = GammaBC(end,1) - 2 ./ (1 + r(end,1)) .* Gamma_plus_collettore;

Gamma_plus_emettitore = gamma * sum(GammaBC(1,qs>0),2);
GammaBC(1,1) = GammaBC(1,1) - 2 ./ (1 + r(1,1)) .* Gamma_plus_emettitore;

GammaBC = GammaBC .* [-1; 1];
% The multiplication for [-1; 1] is to go back to the convention of the
% code, that is, fluxes from left to right

end
