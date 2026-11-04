function Teff = ComputeTeffIon(Ec, mu_c, T, m_ion, M_n, M_b, IonTemp)
%COMPUTETEFFION Effective temperature of ion-neutral collisions, eq. (19) of
%Sieck, Herron, Green, Plasma Chem. Plasma Process. 20, 235 (2000).
%   Ec: electric field at the nodes (V/m), mu_c: ion mobility at the nodes
%   (m^2/Vs), T: gas temperature (K), m_ion, M_n, M_b: masses (amu) of the
%   ion, of the reactive neutral and of the bulk gas. IonTemp = -1 returns T.

if isequal(IonTemp, -1)
    Teff = T*ones(size(Ec));
    return
end

amu = 1.66054e-27;
kB = 1.380649e-23;
vd = mu_c .* abs(Ec);
Teff = T + M_n*amu*vd.^2 * (m_ion + M_b)/(m_ion + M_n) / (3*kB);

end
