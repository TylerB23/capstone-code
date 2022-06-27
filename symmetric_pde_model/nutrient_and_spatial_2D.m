%% Parameters

L = 10; % length of ring where nutrients are distributed and cells exist
nutrient_threshold = 2;
T = 10; % time for simulation

xh = 0.05;  % step size on our discretized ring
x = 0:xh:L;

yh = 0.05; % step size in nutrient amount
y = 0:yh:nutrient_threshold;

k = 0.05;  % step size in time
t = 0:k:T;

