function [Mi] = CreateSatoCurrentMatrix(np,ns,qs,e,r_interfaces,E_L_normalized)
%CREATESATOCURRENTMATRIX Summary of this function goes here
%   Detailed explanation goes here
Gamma2J = CreateGamma2J(np,ns,qs,e);
integral_vector = CreateIntegralVector(diff(r_interfaces));
Mi = 2*pi * integral_vector * spdiags(r_interfaces,0,np+1,np+1) * spdiags(E_L_normalized,0,np+1,np+1) * Gamma2J;
end
