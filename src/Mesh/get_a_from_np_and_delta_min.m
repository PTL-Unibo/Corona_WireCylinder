function [a] = get_a_from_np_and_delta_min(np,delta_min_relative)
%GET_A_FROM_NP_AND_DELTA_MIN Summary of this function goes here
%   Detailed explanation goes here
k = (2 - np) / np;
rhs = 2*delta_min_relative - 1;
fun = @(a) tanh(k*a) - rhs*tanh(a);
a = fzero(fun,[1e-10,10]);
end
