clear; close all; clc;

a = 30;

lambda_num = 200;
lambda = linspace(0.42, 0.68, lambda_num);   % in um
m = water_refractive_index(lambda);

dq = 0.01;
theta1 = 134:dq:148;
intensity1 = water_drop_scattering(a, lambda, theta1, 'SunSize', 0.5);

dq = 0.005;
theta2 = 134:dq:148;
intensity2 = water_drop_scattering(a, lambda, theta2, 'SunSize', 0.5);

%%
figure(1); clf;
subplot(1,2,1);
imagesc(lambda, theta1, intensity1);
subplot(1,2,2);
imagesc(lambda, theta2, intensity2);
