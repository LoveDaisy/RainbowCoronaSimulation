clear; close all; clc;

dq = 0.01;
theta = (137:dq:145)';
a_store = [1000, 700, 450, 300];

%%
figure(1); clf;
set(gcf, 'Position', [50, 500, 900, 450]);

for ai = 1:length(a_store)
    a = a_store(ai);

    lambda_num = max(floor(10000 / a), 50);
    lambda = linspace(0.4, 0.68, lambda_num);   % in um
    dw = lambda(2) - lambda(1);
    sun_spec = colorvis.black_body_radiance(lambda * 1000, 5700);
    sun_spec = sun_spec / sum(sun_spec * dw);

    intensity = water_drop_scattering(a, lambda, theta);
    intensity = intensity ./ sum(intensity * dw * dq) .* sun_spec;
    intensity = intensity ./ prctile(intensity(:), 95) * 5e-3;

    colors = spec_to_rgb([lambda(:) * 1000, intensity'], 'Y', -1);

    axes('Position', [.05, .1+(.79/4+.02)*(ai-1), .92, .79/4]);
    imagesc(theta, 1, reshape(colors, [1, length(theta), 3]));
    set(gca, 'ytick', 1, 'yticklabel', {sprintf('%d', a)}, 'yticklabelrotation', 90, ...
        'tickdir', 'out', 'ticklength', [.006, .015], 'fontsize', 12);
    if ai == 1
        xlabel('Scattering angle (degree)', 'fontsize', 16);
    else
        set(gca, 'xtick', []);
    end
    box off;
    drawnow;
end
axes('Position', [.028, 0.1, 0.0, 0.85], 'Color', 'none');
set(gca, 'ytick', [], 'ycolor', 'none');
ylabel('Radius (um)', 'fontsize', 16, 'color', 'k', 'Visible', 'on');

saveas(gcf, '../img/rainbow_radius_angle.png');
