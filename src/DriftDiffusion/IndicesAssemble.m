function [I_D,J_D,I_Mu,J_Mu,li] = IndicesAssemble(np,ns)
%INDICESASSEMBLED Summary of this function goes here
%   Detailed explanation goes here
I_plus = (2:np)' + (0:(np+1):(np+1)*(ns-1));
J_plus = I_plus - (1:ns);
I_plus = reshape(I_plus,[],1);
J_plus = reshape(J_plus,[],1);
I_minus = I_plus;
J_minus = J_plus + 1;
I_D = [I_plus; I_minus];
J_D = [J_plus; J_minus];
I_Mu = I_plus;
J_Mu = J_plus;

indices = ones((np+1)*ns,1);
excluded_indices = [1;np+1] + (0:(np+1):(np+1)*(ns-1));
indices(reshape(excluded_indices,[],1)) = 0;
li = logical(indices);

end
