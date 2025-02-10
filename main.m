clc;
clear;
close all;

power_gap = [0,150,-200,270,180,150,200,-150,200,270,-180];
time_delay = 0.5; %带二次调频,且每次二次调频响应延迟
%time_delay = 6; %意味着没有二次调频
disp(power_gap)
[frequency, time] = simulate_secondary_frequency(power_gap, time_delay);
disp(frequency)
disp(size(frequency))

%0.008
%0.003
%0.8