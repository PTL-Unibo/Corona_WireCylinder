function [Kelet,ExtraRho,Ke,Add_e,E2,D] = Electrostatics(r_interfaces)

r_nodes = 0.5 * (r_interfaces(1:end-1) + r_interfaces(2:end));
delta = [r_nodes(1)-r_interfaces(1);r_nodes(2:end)-r_nodes(1:end-1);r_interfaces(end)-r_nodes(end)];
coeff = r_interfaces ./ delta;
Vol = (r_interfaces(2:end).^2 - r_interfaces(1:end-1).^2) * 0.5;

np = numel(r_nodes);
i_diag_sx = 2:np; 
j_diag_sx = 1:np-1;
v_diag_sx = -coeff(2:end-1);

i_diag_dx = j_diag_sx;
j_diag_dx = i_diag_sx;
v_diag_dx = v_diag_sx;

i_centro = 2:np-1;
j_centro = i_centro;
v_centro = coeff(3:end-1) + coeff(2:end-2);

i_bcw = 1;
j_bcw = 1;
v_bcw = coeff(2)+coeff(1);

i_bce = np;
j_bce = np;
v_bce = coeff(end)+coeff(end-1);

i = [i_bcw;i_diag_sx';i_centro';i_diag_dx';i_bce];
j = [j_bcw;j_diag_sx';j_centro';j_diag_dx';j_bce];
v = [v_bcw;v_diag_sx;v_centro;v_diag_dx;v_bce];

Kelet = sparse(i,j,v,np,np);
ExtraRho = zeros(np,2);
ExtraRho(1,1) = coeff(1);
ExtraRho(end,end) = coeff(end);

Kelet = Kelet./Vol;
ExtraRho = ExtraRho./Vol;

coeff2 = 1./(r_nodes(2:end) - r_nodes(1:end-1));
coeff2 = [1/(r_nodes(1) - r_interfaces(1)); coeff2; 1/(r_interfaces(end) - r_nodes(end))];

i_sx = 2:np+1;
j_sx = 1:np;
v_sx = coeff2(2:end);

i_dx = 1:np;
j_dx = i_dx;
v_dx = -coeff2(1:end-1);

i = [i_sx';i_dx'];
j = [j_sx';j_dx'];
v = [v_sx;v_dx];

Ke = sparse(i,j,v,np+1,np);
Add_e = zeros(np+1,2);
Add_e(1,1) = coeff2(1);
Add_e(end,end) = -coeff2(end);

E1 = Ke/Kelet;
E2 = E1 * ExtraRho + Add_e;

%-------------------------------------------
row_scale = 1 ./ max(abs(Kelet), [], 2);
D = spdiags(row_scale, 0, size(Kelet,1), size(Kelet,2));

Kelet = D * Kelet;
ExtraRho = D * ExtraRho;
end