function [Gamma2GammaClose] = CreateGamma2GammaClose(np,ns,Vol,r_interfaces)
%CREATEGAMMA2GAMMACLOSE Summary of this function goes here
%   Detailed explanation goes here
I_minus = (1:np)' + (0:np:np*(ns-1));

J_minus = I_minus + (0:ns-1);
I_minus = reshape(I_minus,[],1);
J_minus = reshape(J_minus,[],1);
I_plus = I_minus;
J_plus = J_minus + 1;

S_vol_only = repmat(1./Vol,ns,1);
S_plus = S_vol_only .* repmat(r_interfaces(2:end),ns,1);
S_minus = -S_vol_only .* repmat(r_interfaces(1:end-1),ns,1);

Gamma2GammaClose = sparse([I_minus;I_plus],[J_minus;J_plus],[S_minus;S_plus],np*ns,(np+1)*ns);

end
