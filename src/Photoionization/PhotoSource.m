function [Sph] = PhotoSource(E,s,r_nodes,r_interfaces,Vol,N,Qf,Kph_1,Kph_2,Kph_3,A1,A2,A3,val)
% J. Appl. Phys. 130, 121101 (2021); doi: 10.1063/5.0057856

%---------------------preconditioning-------------------------
row_scale1 = 1 ./ max(abs(Kph_1), [], 2);
row_scale2 = 1 ./ max(abs(Kph_2), [], 2);
row_scale3 = 1 ./ max(abs(Kph_3), [], 2);

Pph_1 = spdiags(row_scale1, 0, size(Kph_1,1), size(Kph_1,2));
Pph_2 = spdiags(row_scale2, 0, size(Kph_2,1), size(Kph_2,2));
Pph_3 = spdiags(row_scale3, 0, size(Kph_3,1), size(Kph_3,2));

Kph_1 = Pph_1 * Kph_1;
Kph_2 = Pph_2 * Kph_2;
Kph_3 = Pph_3 * Kph_3;

%-------------------------------------------------------------
Ec = interp1(r_interfaces,E,r_nodes);
I = (0.03 + (15.7)./(abs(Ec)/N*1e21)) .*Qf .*s;

Sph_1 = Kph_1 \ (Pph_1*(-I .* A1 .* Vol));
Sph_2 = Kph_2 \ (Pph_2*(-I .* A2 .* Vol));
Sph_3 = Kph_3 \ (Pph_3*(-I .* A3 .* Vol));

if val == 1
    Sph = Sph_1 + Sph_2 + Sph_3;
elseif val == -1
    Sph = 0;
end