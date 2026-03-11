function [r_interfaces, delta_max_over_min] = CreateTanhMesh(np,r0,R,a,b)
%CREATETANHMESH Summary of this function goes here
%   Detailed explanation goes here

% "a" can be any positive number
% "b" can be any number smaller than "a" in magnitude

% "a" very small and "b" = 0 corresponds to uniform mesh

tau = linspace(-a-b,a-b,np+1)';
t = tanh(tau);
eta = (t - tanh(-a-b)) / (tanh(a-b) - tanh(-a-b));
r_interfaces = (1 - eta) * r0 + eta * R;

% compute delta_max_over_min to check it's value
delta = diff(eta);
delta_max_over_min = max(delta) / min(delta);

end
