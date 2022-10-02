clear; close all; clc;

a = 50;

lambda_num = 80;
lambda1 = linspace(0.42, 0.68, lambda_num);   % in um
dq = 0.01;
theta1 = 134:dq:148;
tic;
intensity1 = water_drop_scattering(a, lambda1, theta1, 'SunSize', -1);
toc;

lambda_num = 800;
lambda2 = linspace(0.42, 0.68, lambda_num);   % in um
dq = 0.005;
theta2 = 134:dq:148;
tic;
intensity2 = water_drop_scattering(a, lambda2, theta2, 'SunSize', -1);
toc;

%%
figure(1); clf;
subplot(1,2,1);
imagesc(lambda1, theta1, intensity1);
subplot(1,2,2);
imagesc(lambda2, theta2, intensity2);
