function [omega,r,s] = ComputeOmegaTownsend(n_single_column, kr, E, fComputePhotoSource, Const_Omega)
%COMPUTEOMEGATOWNSEND Summary of this function goes here
%   Detailed explanation goes here
persistent Sph

n = reshape(n_single_column,[],3);
omega = zeros(size(n));

n_e   = n(:,1);
n_p  = n(:,2);
n_n  = n(:,3);

r(:,1) = kr(:,1) .* n_e; % ionization
r(:,2) = kr(:,2) .* n_e; % attachment
r(:,3) = kr(:,3) .* n_e .* n_p; % recombination e - p 
r(:,4) = kr(:,3) .* n_n .* n_p; % recombination n - p

s = r(:,1);
if isempty(Sph)
    Sph = fComputePhotoSource(E,s);
end

omega(:,1) = r(:,1) - r(:,2) - r(:,3) + Sph + Const_Omega;
omega(:,2) = r(:,1) - r(:,3) - r(:,4) + Sph + Const_Omega;
omega(:,3) = r(:,2)- r(:,4);

end