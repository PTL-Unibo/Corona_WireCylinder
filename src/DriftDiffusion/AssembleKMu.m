function [KM,u] = AssembleKMu(E,mu,I_Mu,J_Mu,li,qs_col,np,ns)
%ASSEMBLEKMU Summary of this function goes here
%   Detailed explanation goes here
u = qs_col .* repmat(E,ns,1) .* mu;
upwind_indices = u(li) < 0;
KM = sparse(I_Mu,J_Mu+upwind_indices,u(li),(np+1)*ns,np*ns);
end
