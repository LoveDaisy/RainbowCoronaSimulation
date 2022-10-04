clear; close all; clc;

dq = 0.01;
theta_lim = [137, 145];
theta = (theta_lim(1)-0.5:dq:theta_lim(2)+0.5)';
a_store = 10.^(1:.002:3);     % drop radii, in um

lee_diagram = zeros(length(theta), length(a_store), 3);
for ai = length(a_store):-1:1
    a = a_store(ai);
    fprintf('Computing a: %.2f, #%d/%d\n', a, ai, length(a_store));

    lambda_num = max(floor(100000 / a), 500);
    lambda = linspace(0.42, 0.68, lambda_num);   % in um
    dw = lambda(2) - lambda(1);
    sun_spec = colorvis.black_body_radiance(lambda * 1000, 5700);
    sun_spec = sun_spec / sum(sun_spec * dw);

    intensity = water_drop_scattering(a, lambda, theta, 'SunSize', 0.5, 'Parallel', true);
    % energy of each wavelength sums up to that portion of sun light
    intensity = intensity ./ sum(intensity * dw * dq) .* sun_spec;
    intensity = intensity ./ prctile(intensity(:), 95) * 5e-3;

    colors = spec_to_rgb([lambda(:) * 1000, intensity'], 'Y', -1);
    lee_diagram(:, ai, :) = reshape(colors, [], 1, 3);

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
    drawnow;
end

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

postfix = sprintf('a%04d-%04d_q%03d-%03d', min(a_store), max(a_store), theta_lim(1), theta_lim(2));
theta_idx = theta >= theta_lim(1) & theta <= theta_lim(2);
saveas(gcf, sprintf('../out/lee_diagram_%s.png', postfix));
imwrite(uint8(lee_diagram(theta_idx, :, :) * 256), sprintf('../out/lee_diagram_data_%s.png', postfix));