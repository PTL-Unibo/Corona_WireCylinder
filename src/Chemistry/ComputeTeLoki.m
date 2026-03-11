function [Te] = ComputeTeLoki(E, Evals, Tevals, N, val)
%COMPUTETELOKI Summary of this function goes here
%   Detailed explanation goes here
mod_E_Td = abs(E/N*1e21);

if ischar(val) || isstring(val)
    Te = interp1(Evals, Tevals, mod_E_Td);
    Te(mod_E_Td<=Evals(1)) = Tevals(1);
    Te(mod_E_Td>=Evals(end)) = Tevals(end);
else
    Te = ones(size(E)) * val;
end

end
