function [kr, mu, D] = ComputeCoefficientsTownsend(E,Te,T,N,r_interfaces,r_nodes)
%COMPUTECOEFFICIENTTOWNSEND Summary of this function goes here
%   Detailed explanation goes here
Ec = interp1(r_interfaces,E,r_nodes);
Tec = interp1(r_interfaces,Te,r_nodes);

LUT = load("DATA/Loki_E_mu_D_alpha_eta_800Td.csv");

alpha = MyInterp1(LUT(:,1),LUT(:,4),abs(Ec)/N*1e21);
eta = MyInterp1(LUT(:,1),LUT(:,5),abs(Ec)/N*1e21);

kr = zeros(size(Ec,1), 3);
mu = zeros(size(Ec,1)+1, 3);
D = zeros(size(Ec,1)+1, 3);

mu(:,1) = MyInterp1(LUT(:,1),LUT(:,2),abs(E)/N*1e21);
mu(:,2)= (1/N)* min(0.84e23 * (T).^(-0.5), 2.35e12 * (abs(E)/N).^(-0.5));
mu(:,3) = (1/N) * min(0.97e23 * (T).^(-0.5), 3.56e19 * (abs(E)/N).^(-0.1));

mu_ec = interp1(r_interfaces,mu(:,1),r_nodes);

D(:,1) = MyInterp1(LUT(:,1),LUT(:,3),abs(E)/N*1e21);
D(:,[2,3]) = mu(:,[2,3]) .* T / 11600;

kr(:,1) = alpha .* mu_ec .* abs(Ec);
kr(:,2) = eta .* mu_ec .* abs(Ec);
kr(:,3) = 2e-7 * 1e-6;

mu(:,1) = mu(:,1)./N;
D(:,1) = D(:,1)./N;

end