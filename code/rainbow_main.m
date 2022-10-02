clear; close all; clc;

dq = 0.01;
theta = (155:-dq:100)';

img_size = 4096;
[xx, yy] = meshgrid(0:img_size-1, 0:img_size-1);
rr = sqrt((xx - (img_size - 1) / 2).^2 + (yy - (img_size - 1) / 2).^2) / (img_size / 2 * sqrt(2));
rr = rr * 90;
clear xx yy
    
intensity_shaping = exp(min(180 - theta - 45, 0) / 4);
for a = [10, 20, 50, 100, 200, 500, 1000]    % drop radii, in um
    if a <= 30
        lambda_num = 200;
    elseif a <= 80
        lambda_num = 100;
    else
        lambda_num = 50;
    end
    lambda = linspace(0.4, 0.68, lambda_num);   % in um
    sun_spec = colorvis.black_body_radiance(lambda * 1000, 5700);
    sun_spec = sun_spec / sum(sun_spec);

    intensity = water_drop_scattering(a, lambda, theta);
    intensity = intensity ./ max(intensity);

    colors = spec_to_rgb([lambda(:) * 1000, intensity' .* sun_spec' .* intensity_shaping' * 0.02 * lambda_num], 'Y', -1);

    figure(1); clf;
    set(gcf, 'Position', [50, 500, 1200, 100]);
    set(gca, 'Position', [.025, .2, .96, .75]);
    imagesc(theta, 1, reshape(colors, [1, length(theta), 3]));
    saveas(gcf, sprintf('../img/rainbow_mie_%04d.png', a));

    img = interp1(180 - theta, colors, rr, 'linear');
    img(repmat(rr < 0.25, [1, 1, 3])) = 1;

    figure(2); clf;
    imshow(img);
    imwrite(uint8(img * 256), sprintf('../img/rainbow_%04d.png', a));
end
