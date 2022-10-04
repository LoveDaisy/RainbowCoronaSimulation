function lee_diagram_data = generate_lee_diagram(drop_radius, theta, varargin)
p = inputParser;
p.addRequired('drop_radius', @(x) isnumeric(x) && isvector(x));
p.addRequired('theta', @(x) isnumeric(x) && isvector(x));
p.addParameter('SunSize', 0.5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('LambdaRelErr', 0.1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('IntensityFactor', [95, 1], @(x) isnumeric(x) && isvector(x) && length(x) == 2);
p.addParameter('IterateDir', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Debug', false, @(x) islogical(x) && isscalar(x));
p.parse(drop_radius, theta, varargin{:});

dq = abs(theta(2) - theta(1));
theta_lim = [min(theta), max(theta)];
a_num = length(drop_radius);
lee_diagram_data = zeros(length(theta), a_num, 3);
if p.Results.IterateDir > 0
    ai_store = 1:a_num;
else
    ai_store = a_num:-1:1;
end
for ai = ai_store
    a = drop_radius(ai);
    fprintf('Computing a: %.2f, #%d/%d\n', a, ai, a_num);

    [intensity, lambda] = water_drop_scattering(a, [0.42, 0.68], theta, 'SunSize', p.Results.SunSize, ...
        'LambdaRelErr', p.Results.LambdaRelErr, 'AdaptiveLambda', true, 'Parallel', true);
    fprintf('  sample lambda: %d\n', length(lambda));

    dw = lambda(2) - lambda(1);
    sun_spec = colorvis.black_body_radiance(lambda * 1000, 5700);
    sun_spec = sun_spec / sum(sun_spec * dw);

    % energy of each wavelength sums up to that portion of sun light
    intensity = intensity ./ sum(intensity * dw * dq) .* sun_spec;
    intensity = intensity ./ prctile(intensity(:), p.Results.IntensityFactor(1)) * 5e-3 * p.Results.IntensityFactor(2);

    colors = spec_to_rgb([lambda(:) * 1000, intensity'], 'Y', -1);
    lee_diagram_data(:, ai, :) = reshape(colors, [], 1, 3);

    if p.Results.Debug
        axes('Position', [.13, .08, .84, .9]);
        imagesc(log10(drop_radius), theta, lee_diagram_data);
        set(gca, 'ydir', 'reverse', 'ylim', theta_lim);
        axis off;
        axes('Position', [.13, .08, .84, .9], 'Color', 'none');
        set(gca, 'ydir', 'reverse', 'ylim', theta_lim, ...
            'xlim', [min(drop_radius), max(drop_radius)]/1000, 'xscale', 'log', 'tickdir', 'out', ...
            'fontsize', 12);
        box on;
        xlabel('Radius (mm)', 'fontsize', 16);
        ylabel('Scattering angle (degree)', 'fontsize', 16);
        drawnow;
    end
end
end