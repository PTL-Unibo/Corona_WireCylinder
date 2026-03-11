function [indices_boundaries] = CreateIndicesBoundaries(np,ns)
%CREATEINDICESBOUNDARIES Summary of this function goes here
%   Detailed explanation goes here

indices_boundaries_interfaces = [1;np+1] + (0:(np+1):(np+1)*(ns-1));
indices_boundaries_interfaces = reshape(indices_boundaries_interfaces,[],1);

indices_boundaries_nodes = [1;np] + (0:np:np*(ns-1));
indices_boundaries_nodes = reshape(indices_boundaries_nodes,[],1);

indices_boundaries.interfaces = indices_boundaries_interfaces;
indices_boundaries.nodes = indices_boundaries_nodes;

indices_boundaries.interp_left = (indices_boundaries_nodes(1:2:end) + [0,1,2])';
indices_boundaries.interp_right = (indices_boundaries_nodes(2:2:end) + [0,-1,-2])';

end
