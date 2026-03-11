function [integral_vector] = CreateIntegralVector(delta_interfaces)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
integral_vector = zeros(1,numel(delta_interfaces)+1);
integral_vector(1:end-1) = delta_interfaces';
integral_vector(2:end) = integral_vector(2:end) + delta_interfaces';
integral_vector = 0.5 * integral_vector;
end
