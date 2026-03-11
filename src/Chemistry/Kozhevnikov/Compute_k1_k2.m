function [k1,k2] = Compute_k1_k2(E, E_k1_k2, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% converting electric field in Townsend
mod_E_Td = abs(E/N*1e21);

% getting logical indices of where electric field is higher or lower than
% the one of the look-up-table
li_inf = mod_E_Td <= E_k1_k2(1,1); % E(1)
li_SUP = mod_E_Td >= E_k1_k2(end,1); % E(end)

% interpolating look-up-table
k1 = interp1(E_k1_k2(:,1),E_k1_k2(:,2),mod_E_Td); % interp1(E_table,k1_table,E_Td)
k2 = interp1(E_k1_k2(:,1),E_k1_k2(:,3),mod_E_Td); % interp1(E_table,k2_table,E_Td)

% setting k1 values where there is a low electric field to lower value in k1
k1(li_inf) = E_k1_k2(1,2);
k1(li_SUP) = E_k1_k2(end,2);

% setting k2 values where there is a low electric field to lower value in k2
k2(li_inf) = E_k1_k2(1,3);
k2(li_SUP) = E_k1_k2(end,3);

end
