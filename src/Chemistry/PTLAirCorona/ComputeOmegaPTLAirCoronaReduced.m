function [omega,r,s] = ComputeOmegaPTLAirCoronaReduced(n_single_column, kr, E, fComputePhotoSource, Const_Omega)
%COMPUTEOMEGAPTLAIRCORONAREDUCED Compute the kinetic source terms of the reduced PTL-AirCorona
%kinetic scheme (13 reactions). The species vector is the same as in the full scheme
%(e, N2, O2, N2+, O2+, O2-, O, N, O4+); O and N are not involved in any reaction of the reduced
%scheme and keep their initial density. Reaction numbers of the full scheme in brackets.
persistent Sph

n = reshape(n_single_column,[],9);
omega = zeros(size(n));

n_e = n(:,1);
n_N2 = n(:,2);
n_O2 = n(:,3);
n_N2p = n(:,4);
n_O2p = n(:,5);
n_O2m = n(:,6);
n_O4p = n(:,9);

r = zeros(size(n,1), 13);
r(:,1)  = kr(:,1) .* n_e .* n_N2;            % e + N2 -> N2+ + 2e (R1)
r(:,2)  = kr(:,2) .* n_e .* n_O2;            % e + O2 -> O2+ + 2e (R2)
r(:,3)  = kr(:,3) .* n_e .* n_O2 .* n_O2;    % e + 2O2 -> O2- + O2 (R11)
r(:,4)  = kr(:,4) .* n_e .* n_O2 .* n_N2;    % e + O2 + N2 -> O2- + N2 (R12)
r(:,5)  = kr(:,5) .* n_O2m .* n_O2;          % O2- + O2 -> e + 2O2 (R13)
r(:,6)  = kr(:,6) .* n_O2 .* n_N2 .* n_O2p;  % O2+ + O2 + N2 -> O4+ + N2 (R17)
r(:,7)  = kr(:,7) .* n_O2 .* n_O2 .* n_O2p;  % O2+ + 2O2 -> O4+ + O2 (R18)
r(:,8)  = kr(:,8) .* n_O2 .* n_O4p;          % O4+ + O2 -> O2+ + 2O2 (R19)
r(:,9)  = kr(:,9) .* n_O4p .* n_O2m;         % O2- + O4+ -> 3O2 (R22)
r(:,10) = kr(:,10) .* n_O2m .* n_O4p .* n_N2; % O2- + O4+ + N2 -> N2 + 3O2 (R23)
r(:,11) = kr(:,11) .* n_O2m .* n_O4p .* n_O2; % O2- + O4+ + O2 -> 4O2 (R24)
r(:,12) = kr(:,12) .* n_O2m .* n_N2;         % O2- + N2 -> e + O2 + N2 (R25)
r(:,13) = kr(:,13) .* n_N2p .* n_O2;         % N2+ + O2 -> O2+ + N2 (R26)

s = r(:,1) + r(:,2);
if isempty(Sph)
    Sph = fComputePhotoSource(E,s);
end

omega(:,1) = r(:,1) + r(:,2) - r(:,3) - r(:,4) + r(:,5) + r(:,12) + Sph + Const_Omega;
omega(:,2) = -r(:,1) + r(:,13);
omega(:,3) = -r(:,2) - r(:,3) - r(:,4) + r(:,5) - r(:,6) - r(:,7) + r(:,8) + 3*r(:,9) + 3*r(:,10) + 3*r(:,11) + r(:,12) - r(:,13);
omega(:,4) = r(:,1) - r(:,13) + Sph*0.8 + Const_Omega*0.5;
omega(:,5) = r(:,2) - r(:,6) - r(:,7) + r(:,8) + r(:,13) + Sph*0.2 + Const_Omega*0.5;
omega(:,6) = r(:,3) + r(:,4) - r(:,5) - r(:,9) - r(:,10) - r(:,11) - r(:,12) + 1e5;
omega(:,9) = r(:,6) + r(:,7) - r(:,8) - r(:,9) - r(:,10) - r(:,11);

end
