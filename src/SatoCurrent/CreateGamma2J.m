function [Gamma2J] = CreateGamma2J(np,ns,qs,e)
%CREATEGAMMA2J Summary of this function goes here
%   Detailed explanation goes here
I = (1:np+1)';
J = I + (0:(np+1):(np+1)*(ns-1));
S = ones(np+1,ns) .* qs;
Gamma2J = e * sparse(repmat(I,ns,1), reshape(J,[],1), reshape(S,[],1), np+1, (np+1)*ns);
end
