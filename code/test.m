clear; close all; clc;

a = 100;
lambda = 0.65;   % in um
m = 1.332;
dq = 0.01;
% theta = 137:.01:145;
theta = 0:dq:20;

[intensity, Q_sct, Q_ext] = mie_theory_scattering(a, m, lambda, theta);
smoothing_kernel = sqrt(1 - ((-0.25:dq:0.25) / 0.25).^2);
smoothing_kernel = smoothing_kernel / sum(smoothing_kernel);

%%
figure(1); clf;
hold on;
plot(theta, intensity);
plot(theta, conv(intensity, smoothing_kernel, 'same'));
box on;
set(gca, 'yscale', 'log');
