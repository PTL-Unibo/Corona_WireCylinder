function [KD] = AssembleKD(D,delta,I,J,li,np,ns)
%ASSEMBLEKD Summary of this function goes here
%   Detailed explanation goes here
S_plus = D(li) ./ repmat(delta(2:end-1),ns,1);
KD = sparse(I,J,[S_plus;-S_plus],(np+1)*ns,np*ns);
end
