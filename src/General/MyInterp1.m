function [varargout] = MyInterp1(x, V, xq)
% MyInterp1(x, V, xq)
% Automatically uses max and min values of V
% for xq that are out of bounds
%
% - x and xq must be vectors
% - V can be a matrix, with several columns, and
%   the result will be the interpolation on each
%   of the columns
%
% x = [1,2,3,4]
%
% V = [ 2, 7
%       4, 10
%       6, 13
%       8, 16 ]
%
% xq = [0.5, 1.5, 4.5] ------> [ 2, 7
%                                3, 8.5
%                                8, 16  ]
Vq = interp1(x,V,xq);
Vq(xq<x(1),:) = V(1,:) .* ones(sum(xq<x(1)),1);
Vq(xq>x(end),:) = V(end,:) .* ones(sum(xq>x(end)),1);
for k = 1:size(Vq,2)
    varargout{k} = Vq(:,k); %#ok<AGROW>
end
end