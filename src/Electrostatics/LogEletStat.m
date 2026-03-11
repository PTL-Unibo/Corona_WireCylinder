function [E1,E2,Kelet,ExtraRho,K_e,Add_e,D] = LogEletStat(r_interfaces)
%CREATELOGELETSTAT Summary of this function goes here
%   Detailed explanation goes here
r_nodes = 0.5 * (r_interfaces(1:end-1) + r_interfaces(2:end));
np = numel(r_nodes);
Vol = 0.5 * (r_interfaces(2:end).^2 - r_interfaces(1:end-1).^2);

coeff = 1 ./ log(r_nodes(2:end) ./ r_nodes(1:end-1));

Kelet = zeros(np,np);
for i = 1:np-1
    Kelet(i:i+1,i:i+1) = Kelet(i:i+1,i:i+1) + [1,-1;-1,1] * coeff(i); 
end
Kelet(1,1) = Kelet(1,1) +  1/log(r_nodes(1)/r_interfaces(1));
Kelet(end,end) = Kelet(end,end) + 1/log(r_interfaces(end)/r_nodes(end));

ExtraRho = zeros(np,2);
ExtraRho(1,1) = 1/log(r_nodes(1)/r_interfaces(1));
ExtraRho(end,end) = 1/log(r_interfaces(end)/r_nodes(end));

Kelet = Kelet./Vol;
ExtraRho = ExtraRho./Vol;

r_full = [r_interfaces(1); r_nodes; r_interfaces(end)];
coeff = 1 ./ log(r_full(2:end)./r_full(1:end-1)) ./ r_interfaces;

K_e = zeros(np+1,np);
for i = 1:np
    K_e(i:i+1,i) = [-coeff(i); coeff(i+1)];
end

Add_e = zeros(np+1,2);
Add_e(1,1) = coeff(1);
Add_e(end,end) = -coeff(end);

E1 = K_e/Kelet;
E2 = E1 * ExtraRho + Add_e;

Kelet = sparse(Kelet);
ExtraRho = sparse(ExtraRho);
K_e = sparse(K_e);
Add_e = sparse(Add_e);

%------------------------------
row_scale = 1 ./ max(abs(Kelet), [], 2);
D = spdiags(row_scale, 0, size(Kelet,1), size(Kelet,2));
Kelet = D * Kelet;
ExtraRho = D * ExtraRho;

end
