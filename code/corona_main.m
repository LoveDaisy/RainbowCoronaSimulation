clear; close all; clc;

dq = 0.005;
theta = (dq:dq:20)';

img_size = 4096;
[xx, yy] = meshgrid(0:img_size-1, 0:img_size-1);
rr = sqrt((xx - (img_size - 1) / 2).^2 + (yy - (img_size - 1) / 2).^2) / (img_size / 2 * sqrt(2));
rr = rr * 20;
clear xx yy

for a = [5, 6, 7, 8, 10, 15, 20, 30, 40, 50]     % drop radii, in um
    fprintf('Computing drop of %.1fum\n', a);

    lambda_num = max(floor(6000 / a), 30);
    lambda = linspace(0.4, 0.68, lambda_num);   % in um
    dw = lambda(2) - lambda(1);
    sun_spec = colorvis.black_body_radiance(lambda * 1000, 6200);
    sun_spec = sun_spec / sum(sun_spec);

    intensity = water_drop_scattering(a, lambda, theta);
    intensity = intensity ./ sum(intensity * dw * dq) .* sun_spec;
    intensity = intensity ./ prctile(intensity(:), 100) * 2e-3;

    for ev = [0, 2, 4]
        colors = spec_to_rgb([lambda(:) * 1000, intensity' * 2^ev], 'Y', -1);

        figure(1); clf;
        set(gcf, 'Position', [50, 500, 1200, 100]);
        set(gca, 'Position', [.025, .2, .96, .75]);
        imagesc(theta, 1, reshape(colors, [1, length(theta), 3]));
        drawnow;
        saveas(gcf, sprintf('../out/corona_mie_r%04d_e%+2d.png', a, ev));

        img = interp1(theta, colors, rr, 'linear');
        img(repmat(rr < 0.25, [1, 1, 3])) = 1;

        figure(2); clf;
        imshow(img);
        drawnow;
        imwrite(uint8(img * 256), sprintf('../out/corona_r%04d_e%+2d.png', a, ev));
    end
end
