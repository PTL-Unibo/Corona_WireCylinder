function [kr, mu, D] = ComputeCoefficientsHybrid(E,Te,T,N,r_interfaces,r_nodes,MuConst)
%COMPUTECOEFFICIENTSPARENT Compute the reaction coefficients and transport properties according to
%the Parent Model

Ec = interp1(r_interfaces,E,r_nodes);

Tc = T;
Tec = interp1(r_interfaces,Te,r_nodes);

kr = zeros(size(Ec,1), 26);
mu = zeros(size(Ec,1)+1, 9);
D = zeros(size(Ec,1)+1, 9);

kr(:,1) = exp(-0.0105809*log(abs(Ec)/N).^2 - 2.40411e-75*log(abs(Ec)/N).^46) * 1e-6;
kr(:,2) = exp(-0.0102785*log(abs(Ec)/N).^2 - 2.42260e-75*log(abs(Ec)/N).^46) * 1e-6;
kr(:,3) = 2.0e-7*(300./(Tec*11600)).^0.7 * 1e-6;
kr(:,4) = 2.8e-7*(300./(Tec*11600)).^0.5 * 1e-6;
kr(:,5) = 2.0e-7*(300./(Tc)).^0.5 * 1e-6;
kr(:,6) = 2.0e-7*(300./(Tc)).^0.5 * 1e-6;
kr(:,7) = 2.0e-25*(300./Tc).^2.5 * 1e-12;
kr(:,8) = 2.0e-25*(300./Tc).^2.5 * 1e-12;
kr(:,9) = 2.0e-25*(300./Tc).^2.5 * 1e-12;
kr(:,10) = 2.0e-25*(300./Tc).^2.5 * 1e-12;
kr(:,11) = 1.4e-29*(300./(Tec*11600)).*exp(-600./Tc).*exp(700*((Tec*11600)-Tc)./((Tec*11600).*Tc)) * 1e-12;
kr(:,12) = 1.07e-31*(300./(Tec*11600)).^2.*exp(-70./Tc).*exp(1500*((Tec*11600)-Tc)./((Tec*11600).*Tc)) * 1e-12;
kr(:,13) = 8.6e-10*exp(-6030./Tc).*(1-exp(-1570./Tc)) * 1e-6;
kr(:,14) = 1E-11 * 1E-6;
kr(:,15) = 6E-29*(300./(Tec*11600)).^1.5 * 1E-12;
kr(:,16) = 6E-29*(300./(Tec*11600)).^1.5 * 1E-12;
kr(:,17) = 1.4E-6*(300./(Tec*11600)).^0.5 * 1E-6;
kr(:,18) = 1.0E-31 * 1E-12;
kr(:,19) = 2.4E-30 * 1E-12;
kr(:,20) = 2.4E-30 * 1E-12;
kr(:,21) = 1.7E-13 * 1E-6;
kr(:,22) = 1E-7 * (Tec.^(-1.6)) .* exp(-9.8./Tec) * 1E-6;
kr(:,23) = 4.2E-9 * exp(-5.6./Tec) * 1E-6;
kr(:,24) = 1E-7 * 1E-6;
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
