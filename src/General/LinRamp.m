function [v] = LinRamp(t,t_ramp,v0,v_end,t_start)
% t_start is set to 0 as default
if nargin < 5
    t_start = 0;
end
v = v0 + (v_end-v0)*((t-t_start)/t_ramp .* (t>t_start & t<t_start+t_ramp) + (t>=t_start+t_ramp));
end