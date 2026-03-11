function [kr, mu, D] = ComputeCoefficientsParentConst(E,Te,T,N,r_interfaces,r_nodes,MuConst)
%COMPUTECOEFFICIENTSPARENT Compute the reaction coefficients and transport properties according to
%the Parent Model

Ec = interp1(r_interfaces,E,r_nodes);

Tc = T;
Tec = interp1(r_interfaces,Te,r_nodes);

kr = zeros(size(Ec,1), 13);
mu = zeros(size(Ec,1)+1, 4);
D = zeros(size(Ec,1)+1, 4);

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

if isequal(MuConst, -1)
    mu(:,1) = (1/N) * 3.74e19 * exp(33.5 * (log(Te*11600)).^(-0.5));
    mu(:,2) = (1/N) * min(0.75e23 * (T).^(-0.5), 2.03e12 * (abs(E)/N).^(-0.5));
    mu(:,3) = (1/N) * min(1.18e23 * (T).^(-0.5), 3.61e12 * (abs(E)/N).^(-0.5));
    mu(:,4) = (1/N) * min(0.97e23 * (T).^(-0.5), 3.56e19 * (abs(E)/N).^(-0.1));
elseif isvector(MuConst)
    mu(:,1) = MuConst(1);
    mu(:,2) = MuConst(2);
    mu(:,3) = MuConst(3);
    mu(:,4) = MuConst(4);
end

D(:,1) = mu(:,1) .* Te;
D(:,[2,3,4]) = mu(:,[2,3,4]) .* T / 11600;

end
