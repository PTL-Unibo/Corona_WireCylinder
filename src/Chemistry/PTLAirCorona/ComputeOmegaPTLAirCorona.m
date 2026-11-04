function [omega,r,s] = ComputeOmegaPTLAirCorona(n_single_column, kr, E, fComputePhotoSource, Const_Omega)
%COMPUTEOMEGAPTLAIRCORONA Compute the kinetic source terms of the PTL-AirCorona
%kinetic scheme (9 species, 26 reactions)
persistent Sph

n = reshape(n_single_column,[],9);
omega = zeros(size(n));

n_e = n(:,1);
n_N2 = n(:,2);
n_O2 = n(:,3);
n_N2p = n(:,4);
n_O2p = n(:,5);
n_O2m = n(:,6);
n_O = n(:,7);
n_N = n(:,8);
n_O4p = n(:,9);

r(:,1)  = kr(:,1) .* n_e .* n_N2;
r(:,2)  = kr(:,2) .* n_e .* n_O2;
r(:,3)  = kr(:,3) .* n_e .* n_O2p;
r(:,4)  = kr(:,4) .* n_e .* n_N2p;
r(:,5)  = kr(:,5) .* n_O2m .* n_N2p;
r(:,6)  = kr(:,6) .* n_O2m .* n_O2p;
r(:,7)  = kr(:,7) .* n_O2m .* n_N2p .* n_N2;
r(:,8)  = kr(:,8) .* n_O2m .* n_O2p .* n_N2;
r(:,9)  = kr(:,9) .* n_O2m .* n_N2p .* n_O2;
r(:,10) = kr(:,10) .* n_O2m .* n_O2p .* n_O2;
r(:,11) = kr(:,11) .* n_e .* n_O2 .* n_O2;
r(:,12) = kr(:,12) .* n_e .* n_O2 .* n_N2;
r(:,13) = kr(:,13) .* n_O2m .* n_O2;
r(:,14) = kr(:,14) .* n_e .* n_e .* n_O2p;
r(:,15) = kr(:,15) .* n_e .* n_O4p;
r(:,16) = kr(:,16) .* n_e .* n_O2 .* n_O;
r(:,17) = kr(:,17) .* n_O2 .* n_N2 .* n_O2p;
r(:,18) = kr(:,18) .* n_O2 .* n_O2 .* n_O2p;
r(:,19) = kr(:,19) .* n_O2 .* n_O4p;
r(:,20) = kr(:,20) .* n_e .* n_N2;
r(:,21) = kr(:,21) .* n_e .* n_O2;
r(:,22) = kr(:,22) .* n_O4p .* n_O2m;
r(:,23) = kr(:,23) .* n_O2m .* n_O4p .* n_N2;
r(:,24) = kr(:,24) .* n_O2m .* n_O4p .* n_O2;
r(:,25) = kr(:,25) .* n_O2m .* n_N2;
r(:,26) = kr(:,26) .* n_N2p .* n_O2;

s = r(:,1) + r(:,2);
if isempty(Sph)
    Sph = fComputePhotoSource(E,s);
end

omega(:,1) = -r(:,11) + r(:,1) + r(:,2) - r(:,4) - r(:,3) + r(:,13) - r(:,12) - r(:,14) - r(:,15) - r(:,16) + r(:,25) + Sph + Const_Omega;
omega(:,2) = -r(:,1) + r(:,4) + r(:,5) + r(:,7) + r(:,9) - r(:,20) + r(:,26);
omega(:,3) = -r(:,11) - r(:,2) + r(:,3) + r(:,5) + 2*r(:,6) + r(:,13) - r(:,12) + r(:,7) + 2*r(:,8) + r(:,9) + 2*r(:,10) + r(:,14) + 2*r(:,15) - r(:,16) - r(:,17) - r(:,18) + r(:,19) - r(:,21) + 3*r(:,22) + 3*r(:,23) + 3*r(:,24) + r(:,25) - r(:,26);
omega(:,4) = r(:,1) - r(:,4) - r(:,5) - r(:,7) - r(:,9) - r(:,26) + Sph*0.8 + Const_Omega*0.5;
omega(:,5) = r(:,2) - r(:,3) - r(:,6) - r(:,8) - r(:,10) - r(:,14) - r(:,17) - r(:,18) + r(:,19) + r(:,26) + Sph*0.2 + Const_Omega*0.5;
omega(:,6) = r(:,11) - r(:,5) - r(:,6) - r(:,13) + r(:,12) - r(:,7) - r(:,8) - r(:,9) - r(:,10) + r(:,16) - r(:,22) - r(:,23) - r(:,24) - r(:,25) + 1e5;
omega(:,7) = 2*r(:,21);
omega(:,8) = 2*r(:,20);
omega(:,9) = - r(:,15) + r(:,17) + r(:,18) - r(:,19) - r(:,22) - r(:,23) - r(:,24);

end
