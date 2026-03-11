function [Jpattern] = CreateJacobianSparsityPattern(np, ns, qs)
%CREATEJACOBIANSPARSITYPATTERN Summary of this function goes here
%   Detailed explanation goes here
%assemblo la parte delle n
J_N = (1 : np : np*ns)' + (0 : np-1);
J_N = repmat(reshape(J_N, [], 1), ns, 1);
I_N = (1 : np*ns)';
I_N = repelem(I_N(:), ns);
I_N_diagSX = reshape(((2 : np)' + (0 : np : np*(ns-1))), [], 1);
J_N_diagSX = I_N_diagSX - 1;
I_N_diagDX = J_N_diagSX;
J_N_diagDX = I_N_diagSX;
%assemblo i phi diagonali
J_phi_diag = (1 : np*ns)';
I_phi_diag = repmat((np*ns + 1 : np*(ns + 1))', ns, 1);
k = 1 : np : np*ns;
for i = 1 : length(qs)
     if qs(i) == 0
     I_phi_diag(k(i) : np*i) = 0;
     J_phi_diag(k(i) : np*i) = 0;
     end
end
J_phi_diag = J_phi_diag(J_phi_diag > 0);
I_phi_diag = I_phi_diag(I_phi_diag > 0);
%assemblo i phi verticali
J_phi = ns*np + 1 : np*(ns + 1);
J_phiSx = repelem(ns*np + 1, 2*(ns + 1), 1);
J_phiDx = repelem(np*(ns + 1), 2*(ns + 1), 1);
J_phi_Central = reshape(repelem((J_phi(2) : J_phi(end - 1))', 3*(ns + 1)), [], 1);
I_phiSx = 1 : np : np*(ns + 1);
I_phiSx = reshape([I_phiSx; I_phiSx + 1], [], 1);
I_phiDx = np - 1 : np : np*(ns + 1);
I_phiDx = reshape([I_phiDx; I_phiDx + 1], [], 1);
I_phi_Central = 1 : np : np*(ns + 1);
I_phi_Central = reshape([I_phi_Central; I_phi_Central + 1; I_phi_Central + 2], [], 1);
I_phi_Central = reshape(I_phi_Central + repmat((0 : np - 3), 3*(ns + 1), 1), [], 1);
I_tot = [I_N; I_N_diagSX; I_N_diagDX; I_phi_diag; I_phiSx; I_phi_Central; I_phiDx];
J_tot = [J_N; J_N_diagSX; J_N_diagDX; J_phi_diag; J_phiSx; J_phi_Central; J_phiDx];
Jpattern = sparse(I_tot, J_tot, 1, np*(ns + 1), np*(ns + 1));
end
