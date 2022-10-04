clear; close all; clc;

a = 100;

lambda = 0.65;
dq = 0.01;
m = 1.332;
theta = 137:dq:145;

tic;
intensity = mie_scattering(a, m, lambda, theta);
toc;

%%
figure(1); clf;
plot(theta, intensity);
set(gca, 'yscale', 'linear');
