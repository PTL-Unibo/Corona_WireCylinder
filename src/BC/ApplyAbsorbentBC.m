function [GammaBC] = ApplyAbsorbentBC(n_boundary,u_boundary,ns)
%APPLYABSORBENTBC Summary of this function goes here
%   Detailed explanation goes here
n_boundary = reshape(n_boundary,2,ns);
u_boundary = reshape(u_boundary,2,ns);
GammaBC = n_boundary .* u_boundary;
GammaBC = reshape(GammaBC,[],1);
end
