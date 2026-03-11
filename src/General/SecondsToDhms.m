function [d, h, m, s] = SecondsToDhms(time_in_seconds, print_flag)
% SecondsToDhms Converts "time_in_seconds" to days : hours : minutes : seconds
arguments
    time_in_seconds
    print_flag char {mustBeMember(print_flag,{'print','no'})} = "no"
end
s = round(mod(time_in_seconds, 60));
minutes = fix(time_in_seconds/60);
m = mod(minutes, 60);
hours = fix(minutes/60);
h = mod(hours, 24);
d = fix(hours/24);
if print_flag == "print"
    fprintf("%dd : %dh : %dm : %ds\n", d, h, m, s);
end
end

