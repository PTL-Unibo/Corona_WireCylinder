function [omega,r,s] = ComputeOmegaParentConst(n_single_column, kr, E, N, fComputePhotoSource, Const_Omega)
%COMPUTEOMEGAPARENT Summary of this function goes here
%   Detailed explanation goes here
persistent Sph

n = reshape(n_single_column,[],4);
omega = zeros(size(n));

n_e   = n(:,1);
n_N2  = 0.7884*N;
n_O2  = 0.2116*N;
n_N2p = n(:,2);
n_O2p = n(:,3);
n_O2m = n(:,4);

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

s = r(:,1) + r(:,2);
if isempty(Sph)
    Sph = fComputePhotoSource(E,s);
end

omega(:,1) = -r(:,11) + r(:,1) + r(:,2) - r(:,4) - r(:,3) + r(:,13) - r(:,12) + Sph + Const_Omega;
omega(:,2) = r(:,1) - r(:,4) - r(:,5) - r(:,7) - r(:,9) + Sph*0.5 + Const_Omega*0.5;
omega(:,3) = r(:,2) - r(:,3) - r(:,6) - r(:,8) - r(:,10) + Sph*0.5 + Const_Omega*0.5;
omega(:,4) = r(:,11) - r(:,5) - r(:,6) - r(:,13) + r(:,12) - r(:,7) - r(:,8) - r(:,9) - r(:,10) + 1e5; 

end
