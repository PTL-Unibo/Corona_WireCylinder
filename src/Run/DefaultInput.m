function [opts] = DefaultInput(opts)
%DefaultFROG Summary of this function goes here
%   Detailed explanation goes here
opts.a = 0; 
opts.b = 0;
opts.gamma_II = 1e-3;
opts.Const_Omega = 0;
opts.Poisson = 'Lin';
opts.Te = 'Loki';
opts.MuConst = -1;
opts.T_GAS = 300;
opts.IonTemp = -1; % 1: effective ion temperature (Sieck et al. 2000), -1: T_GAS
opts.RelTol = 1e-2; 
opts.OutFcn = 'bar';
opts.mesh = -1;

end