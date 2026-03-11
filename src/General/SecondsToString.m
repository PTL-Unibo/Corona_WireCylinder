function [time_string] = SecondsToString(time_in_seconds)
% SecondsToString Converts "time_in_seconds" to a string indicating the time

[d, h, m, s] = SecondsToDhms(time_in_seconds);
time_string = "";
first_char = "";
zero_flag = true;

if d ~= 0
    time_string = time_string + first_char + num2str(d) + "d";
    first_char = " ";
    zero_flag = false;
end
if h ~= 0
    time_string = time_string + first_char + num2str(h) + "h";
    first_char = " ";
    zero_flag = false;
end
if m ~= 0
    time_string = time_string + first_char + num2str(m) + "m";
    first_char = " ";
    zero_flag = false;
end
if s ~= 0
    time_string = time_string + first_char + num2str(s) + "s";
    first_char = " ";
    zero_flag = false;
end
if zero_flag
    time_string = "0s";
end
 
end

 