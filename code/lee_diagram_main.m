clear; close all; clc;

dq = 0.01;
theta = (137:dq:145)';
a_store = 10.^(1:.002:3);     % drop radii, in um

lee_diagram = zeros(length(theta), length(a_store), 3);
for ai = 1:length(a_store)
    a = a_store(ai);
    fprintf('Computing a: %.2f, #%d/%d\n', a, ai, length(a_store));
    if a < 200
        lambda_num = floor(6000 / a);
    else
        lambda_num = 30;
    end
    lambda = linspace(0.42, 0.68, lambda_num);   % in um
    sun_spec = colorvis.black_body_radiance(lambda * 1000, 5700);
    sun_spec = sun_spec / sum(sun_spec);

    intensity = water_drop_scattering(a, lambda, theta, 'SunSize', 0.5);
%     intensity = intensity ./ max(intensity(:));

    colors = spec_to_rgb([lambda(:) * 1000, intensity' .* sun_spec' * 2 * lambda_num], 'Y', -1);
    lee_diagram(:, ai, :) = reshape(colors, [], 1, 3);

    figure(1); clf;
    set(gcf, 'Position', [500, 20, 500, 800]);
    set(gca, 'Position', [.1, .05, .8, .9]);
    imagesc(log10(a_store), theta, lee_diagram);

    xtick = [1, 2, 3];
    xtick_label = cell(length(xtick), 1);
    for i = 1:length(xtick)
        xtick_label{i} = sprintf('%.1f', 10^xtick(i));
    end
    set(gca, 'xticklabel', xtick_label, 'xtick', xtick, 'ydir', 'reverse');
    drawnow;
end

%%
figure(1); clf;
set(gcf, 'Position', [500, 20, 500, 800]);
set(gca, 'Position', [.1, .05, .8, .9]);
imagesc(log10(a_store), theta, lee_diagram);

xtick = [1, 2, 3];
xtick_label = cell(length(xtick), 1);
for i = 1:length(xtick)
    xtick_label{i} = sprintf('%.1f', 10^xtick(i));
end
set(gca, 'xticklabel', xtick_label, 'xtick', xtick, 'ydir', 'reverse');
axis ij;
drawnow;
