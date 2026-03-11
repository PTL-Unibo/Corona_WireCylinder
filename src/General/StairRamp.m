function [t,v,time_instants] = StairRamp(t_rest,t_ramp,v_array)
v = repelem(v_array,2);
N = numel(v_array);
c1 = repelem(0:(N-1),2);
c2 = [0:(N-1); 1:N]; 
c1 = c1(:);
c2 = c2(:);
t = c1 * t_ramp + c2 * t_rest;
time_instants = 0.999*t_rest + (0:N-1)*(t_rest+t_ramp);
end
