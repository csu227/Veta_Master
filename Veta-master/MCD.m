function [upperbound,lowerbound] = MCD(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
MCD = mean(abs(diff(x)));
avg = mean(x);

upperbound = avg + MCD * 2.66; % used to find offset
lowerbound = avg - MCD * 2.66; %used to find onset
end

