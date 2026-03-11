function [time_string] = GetTimeString(t)
%GETTIMESTRING Summary of this function goes here
%   Detailed explanation goes here

[h,m,s] = hms(t);
[Y,M,D] = ymd(t);
time_string = num2str(Y) + "_" + num2str(M) + "_" + num2str(D) + ...
    "_" + num2str(h) + "_" + num2str(m) + "_" + num2str(round(s));

end
