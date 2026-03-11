function [c1,c2,c3,A1,A2,A3] = PhotoCoefficients(p_O2)

c1 = (0.0533*p_O2*1e2)^2;
c2 = (0.1460*p_O2*1e2)^2;
c3 = (0.4886*p_O2*1e2)^2;

A1 = 1.986e-4*(p_O2^2)*1e4;
A2 = 0.0051*(p_O2^2)*1e4;
A3 = 0.4886*(p_O2^2)*1e4;

end