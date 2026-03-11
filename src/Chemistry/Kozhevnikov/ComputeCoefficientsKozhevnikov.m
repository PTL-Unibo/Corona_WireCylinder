function [kr, mu, D] = ComputeCoefficientsKozhevnikov(E,Te,T,fCompute_k1_k2,N,r_interfaces,r_nodes,MuConst)
    %np,E_look_up_table,k1_look_up_table,k2_look_up_table,N
% ComputeCoefficientsKozhevnikov Summary of this function goes here
%   Detailed explanation goes here

Ec = interp1(r_interfaces,E,r_nodes);
Tc = T;
Tec = interp1(r_interfaces,Te,r_nodes);

kr = zeros(size(Ec,1), 26);
mu = zeros(size(Ec,1)+1, 9);
D = zeros(size(Ec,1)+1, 9);

[k1,k2] = fCompute_k1_k2(Ec);

kr(:,1) = k1;
kr(:,2) = k2;
kr(:,3) = 1E-19*(300./(Tec*11600)).^4.5 * 1e-12;
kr(:,4) = 1E-11 * 1E-6;
kr(:,5) = 6E-29*(300./(Tec*11600)).^1.5 * 1E-12;
kr(:,6) = 6E-29*(300./(Tec*11600)).^1.5 * 1E-12;
kr(:,7) = 6E-29*(300./(Tec*11600)).^1.5 * 1E-12;
kr(:,8) = 6E-29*(300./(Tec*11600)).^1.5 * 1E-12;
kr(:,9) = 1.4E-6*(300./(Tec*11600)).^0.5 * 1E-6;
kr(:,10) = 2E-7*(300./(Tec*11600)).^1.5 * 1E-6;
kr(:,11) = 2.8E-7*(300./(Tec*11600)).^0.5 * 1E-6;
kr(:,12) = 1.9E-30*(300./(Tec*11600)).*exp(7/3*(1-(300./(Tec*11600)))) * 1E-12;
kr(:,13) = 2.2E-18 * 1E-12;
kr(:,14) = 8.5E-32 * (300./(Tec*11600)).^2 .* exp(5*(1-(300./(Tec*11600)))) * 1E-12;
kr(:,15) = 2.2E-19 * 1E-6;
kr(:,16) = 1.0E-31 * 1E-12;
kr(:,17) = 2.4E-30 * 1E-12;
kr(:,18) = 2.4E-30 * 1E-12;
kr(:,19) = 1.7E-13 * 1E-6;
kr(:,20) = 1E-7 * (Tec.^(-1.6)) .* exp(-9.8./Tec) * 1E-6;
kr(:,21) = 4.2E-9 * exp(-5.6./Tec) * 1E-6;
kr(:,22) = 1E-7 * 1E-6;
kr(:,23) = 2E-25 * 1E-12;
kr(:,24) = 2E-25 * 1E-12;
kr(:,25) = 2E-25 * 1E-12;
kr(:,26) = 2E-25 * 1E-12;

if isequal(MuConst, -1)
    mu(:,1) = (1/N) * 3.74e19 * exp(33.5 * (log(Te*11600)).^(-0.5));
    mu(:,4) = (1/N) * min(0.75e23 * (T).^(-0.5), 2.03e12 * (abs(E)/N).^(-0.5));
    mu(:,5) = (1/N) * min(1.18e23 * (T).^(-0.5), 3.61e12 * (abs(E)/N).^(-0.5));
    mu(:,6) = (1/N) * min(0.97e23 * (T).^(-0.5), (3.56e19 * (abs(E)/N).^(-0.1))*1.10);
    mu(:,9) = (1/N) * min(1.18e23 * (T).^(-0.5), 3.61e12 * (abs(E)/N).^(-0.5)) / 1.37; % 1.37 from Sakiyama
elseif isvector(MuConst)
    mu(:,1) = MuConst(1);
    mu(:,4) = MuConst(2);
    mu(:,5) = MuConst(3);
    mu(:,6) = MuConst(4);
    mu(:,9) = MuConst(5);
end

D(:,1) = mu(:,1) .* Te;
D(:,[4,5,6,9]) = mu(:,[4,5,6,9]) .* T / 11600;
D(:,2) = 2.1e-5;
D(:,3) = 2.1e-5;
D(:,7) = 3.2e-5;
D(:,8) = 2.9e-5;

end
