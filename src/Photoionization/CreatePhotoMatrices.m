function [Kph_1,Kph_2,Kph_3,ExtraPH] = CreatePhotoMatrices(np,r_interfaces,delta,Vol,c1,c2,c3)

coeff = r_interfaces ./ delta;

i_c = 1:np;
j_c = i_c;
v_c1 = -coeff(2:end) - coeff(1:end-1) - c1*Vol;
v_c2 = -coeff(2:end) - coeff(1:end-1) - c2*Vol;
v_c3 = -coeff(2:end) - coeff(1:end-1) - c3*Vol;

i_dx = 2:np;
j_dx = 1:np-1;
v_dx = coeff(2:end-1);

i_sx = 1:np-1;
j_sx = 2:np;
v_sx = coeff(2:end-1);

i = [i_c'; i_dx'; i_sx'];
j = [j_c'; j_dx'; j_sx'];
v1 = [v_c1; v_dx; v_sx];
v2 = [v_c2; v_dx; v_sx];
v3 = [v_c3; v_dx; v_sx];

Kph_1 = sparse(i,j,v1,np,np);
Kph_2 = sparse(i,j,v2,np,np);
Kph_3 = sparse(i,j,v3,np,np);

ExtraPH = zeros(np,2);
ExtraPH(1,1) = coeff(1);
ExtraPH(end,end) = coeff(end);

end