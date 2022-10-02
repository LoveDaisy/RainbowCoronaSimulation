clear; close all; clc;

dq = 0.01;
theta = (180:-dq:90)';

img_size = 4096;
[xx, yy] = meshgrid(0:img_size-1, 0:img_size-1);
rr = sqrt((xx - (img_size - 1) / 2).^2 + (yy - (img_size - 1) / 2).^2) / (img_size / 2 * sqrt(2));
rr = rr * 90;
clear xx yy

for a = [10, 20, 50, 100, 200, 500, 1000]    % drop radii, in um
    fprintf('Computing drop of %.1fum\n', a);

    lambda_num = max(floor(6000 / a), 30);
    lambda = linspace(0.4, 0.68, lambda_num);   % in um
    sun_spec = colorvis.black_body_radiance(lambda * 1000, 6200);
    sun_spec = sun_spec / sum(sun_spec);

    intensity = water_drop_scattering(a, lambda, theta);

    colors = spec_to_rgb([lambda(:) * 1000, intensity' .* sun_spec' * 5 * lambda_num], 'Y', -1);

    figure(1); clf;
    set(gcf, 'Position', [50, 500, 1200, 100]);
    set(gca, 'Position', [.025, .2, .96, .75]);
    imagesc(theta, 1, reshape(colors, [1, length(theta), 3]));
    drawnow;
    saveas(gcf, sprintf('../img/rainbow_mie_%04d.png', a));

    img = interp1(180 - theta, colors, rr, 'linear');
    img(repmat(rr < 0.25, [1, 1, 3])) = 1;

    figure(2); clf;
    imshow(img, 'InitialMagnification', 'fit');
    drawnow;
    imwrite(uint8(img * 256), sprintf('../img/rainbow_%04d.png', a));
end
