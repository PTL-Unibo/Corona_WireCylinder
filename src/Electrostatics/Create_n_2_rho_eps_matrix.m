function [N2RhoEps] = Create_n_2_rho_eps_matrix(np,ns,qs,e,eps0)
%CREATE_N_2_RHO_MATRIX Summary of this function goes here
%   Detailed explanation goes here
I = repmat((1:np)',ns,1);
J = reshape((1:np)' + (0:np:np*(ns-1)), [], 1);
S = reshape(ones(np,ns).*qs,[],1) * e / eps0;
N2RhoEps = sparse(I,J,S,np,np*ns);
end
