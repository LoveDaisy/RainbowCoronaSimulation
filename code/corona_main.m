clear; close all; clc;

lambda = linspace(0.4, 0.68, 50);   % in um
sun_spec = colorvis.black_body_radiance(lambda * 1000, 5700);
sun_spec = sun_spec / sum(sun_spec);

dq = 0.005;
theta = (dq:dq:20)';

img_size = 4096;
[xx, yy] = meshgrid(0:img_size-1, 0:img_size-1);
rr = sqrt((xx - (img_size - 1) / 2).^2 + (yy - (img_size - 1) / 2).^2) / (img_size / 2 * sqrt(2));
rr = rr * 20;
clear xx yy

for a = [5, 10, 20, 30, 40, 50]
    fprintf('Corona simulation for a: %.1f\n', a);
    intensity = water_drop_scattering(a, lambda, theta);
    intensity = intensity ./ max(intensity);

    colors = spec_to_rgb([lambda(:) * 1000, intensity' .* sun_spec' * 0.5], 'Y', -1);

    figure(1); clf;
    set(gcf, 'Position', [50, 500, 1200, 100]);
    set(gca, 'Position', [.025, .2, .96, .75]);
    imagesc(theta, 1, reshape(colors, [1, length(theta), 3]));
    saveas(gcf, sprintf('../img/corona_mie_%04d.png', a));

    img = interp1(theta, colors, rr, 'linear');
    img(repmat(rr < 0.25, [1, 1, 3])) = 1;

    figure(2); clf;
    imshow(img);
    imwrite(uint8(img * 256), sprintf('../img/corona_%04d.png', a));
end
