clear; close all; clc;

dq = 0.02;
theta_lim = [170, 180];
theta = (theta_lim(1)-0.5:dq:theta_lim(2)+0.5)';
a_store = 10.^(0:.002:log10(21));     % drop radii, in um

figure(1); clf;
set(gcf, 'Position', [500, 20, 500, 800]);

lee_diagram = generate_lee_diagram(a_store, theta, 'WSampleFactor', 10, 'Debug', true);

%%
figure(1); clf;
set(gcf, 'Position', [500, 20, 500, 800]);
axes('Position', [.13, .08, .84, .9]);
imagesc(log10(a_store), theta, lee_diagram);
set(gca, 'ydir', 'reverse', 'ylim', theta_lim);
axis off;
axes('Position', [.13, .08, .84, .9], 'Color', 'none');
set(gca, 'ydir', 'reverse', 'ylim', theta_lim, ...
    'xlim', [min(a_store), max(a_store)]/1000, 'xscale', 'log', 'tickdir', 'out', ...
    'fontsize', 12);
box on;
xlabel('Radius (mm)', 'fontsize', 16);
ylabel('Scattering angle (degree)', 'fontsize', 16);
drawnow;

postfix = sprintf('a%04d-%04d_q%03d-%03d', round(min(a_store)), round(max(a_store)), theta_lim(1), theta_lim(2));
theta_idx = theta >= theta_lim(1) & theta <= theta_lim(2);
saveas(gcf, sprintf('../out/lee_diagram_%s.png', postfix));
imwrite(uint8(lee_diagram(theta_idx, :, :) * 256), sprintf('../out/lee_diagram_data_%s.png', postfix));