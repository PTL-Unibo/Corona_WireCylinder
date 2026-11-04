function [kr, mu, D] = ComputeCoefficientsPTLAirCoronaReduced(E,Te,T,fCompute_k1_k2,LUT,N,r_interfaces,r_nodes,MuConst)
%COMPUTECOEFFICIENTSPTLAIRCORONAREDUCED Compute the reaction coefficients and transport properties
%of the reduced PTL-AirCorona kinetic scheme (13 reactions, charged species e, N2+, O2+, O2-, O4+).
%   It is obtained from the full PTL-AirCorona scheme by removing reactions R3-R10, R14-R16, R20
%   and R21. Column j of kr corresponds to reaction R(j) of the reduced scheme; the number of the
%   same reaction in the full scheme is reported in the comments. Transport properties are the same
%   as in the full scheme.

Ec = interp1(r_interfaces,E,r_nodes);

Tc = T;
Tec = interp1(r_interfaces,Te,r_nodes);

kr = zeros(size(Ec,1), 13);
mu = zeros(size(Ec,1)+1, 9);
D = zeros(size(Ec,1)+1, 9);

if isequal(MuConst, -1)
    mu(:,1) = MyInterp1(LUT(:,1),LUT(:,2),abs(E)/N*1e21) / N; % electron mobility from the LoKI look-up table
    mu(:,4) = (1/N) * min(0.75e23 * (T).^(-0.5), 2.03e12 * (abs(E)/N).^(-0.5));
    mu(:,5) = (1/N) * min(1.18e23 * (T).^(-0.5), 3.61e12 * (abs(E)/N).^(-0.5));
    mu(:,6) = (1/N) * min(0.97e23 * (T).^(-0.5), 3.56e19 * (abs(E)/N).^(-0.1));
    mu(:,9) = (1/N) * min(1.18e23 * (T).^(-0.5), 3.61e12 * (abs(E)/N).^(-0.5)) / 1.37; % 1.37 from Sakiyama
elseif isvector(MuConst)
    mu(:,1) = MuConst(1);
    mu(:,4) = MuConst(2);
    mu(:,5) = MuConst(3);
    mu(:,6) = MuConst(4);
    mu(:,9) = MuConst(5);
end

if isequal(MuConst, -1)
    D(:,1) = MyInterp1(LUT(:,1),LUT(:,3),abs(E)/N*1e21) / N; % electron diffusion from the LoKI look-up table
else
    D(:,1) = mu(:,1) .* Te;
end
D(:,[4,5,6,9]) = mu(:,[4,5,6,9]) .* T / 11600;
D(:,2) = 2.1e-5;
D(:,3) = 2.1e-5;
D(:,7) = 3.2e-5;
D(:,8) = 2.9e-5;

% ionization rates from the look-up table of the 9-species model
[kr(:,1), kr(:,2)] = fCompute_k1_k2(Ec);                       % R1, R2 (full: R1, R2)
kr(:,3) = 1.4e-29*(300./(Tec*11600)).*exp(-600./Tc).*exp(700*((Tec*11600)-Tc)./((Tec*11600).*Tc)) * 1e-12; % full: R11
kr(:,4) = 1.07e-31*(300./(Tec*11600)).^2.*exp(-70./Tc).*exp(1500*((Tec*11600)-Tc)./((Tec*11600).*Tc)) * 1e-12; % full: R12
kr(:,5) = 1.24e-11*exp(-(179./(8.8+abs(Ec)/N*1e21)).^2) * 1E-6; % full: R13, Pancheshnyi (2013), DOI: 10.1088/0022-3727/46/15/155201
kr(:,6) = 2.4E-30 * 1E-12;                                     % full: R17
kr(:,7) = 2.4E-30 * 1E-12;                                     % full: R18
kr(:,8) = 1.7E-13 * 1E-6;                                      % full: R19
kr(:,9) = 1.0E-7 * 1E-6;                                       % full: R22, Kozhevnikov et al. (2023), Table 1, R19
kr(:,10) = 2.0E-25 * 1E-12;                                    % full: R23, Kozhevnikov et al. (2023), Table 1, R20
kr(:,11) = 2.0E-25 * 1E-12;                                    % full: R24, Kozhevnikov et al. (2023), Table 1, R20
kr(:,12) = 1.24e-11*exp(-(179./(8.8+abs(Ec)/N*1e21)).^2) * 1E-6; % full: R25, Pancheshnyi (2013)
kr(:,13) = 1E-11 * 1E-6;                                       % full: R26, N2+ + O2 -> N2 + O2+, charge transfer (Kozhevnikov)

end
