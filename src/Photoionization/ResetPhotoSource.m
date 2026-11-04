function status = ResetPhotoSource(~,~,~)
%RESETPHOTOSOURCE Clear the persistent photoionization source stored in the
%ComputeOmega* functions, so that it is recomputed at the next call.
%   It can also be used as an ode OutputFcn (it always returns status = 0).
clear ComputeOmegaKozhevnikov
clear ComputeOmegaTownsend
clear ComputeOmegaParent
clear ComputeOmegaParentConst
clear ComputeOmegaPTLAirCorona
clear ComputeOmegaPTLAirCoronaReduced
status = 0;
end
