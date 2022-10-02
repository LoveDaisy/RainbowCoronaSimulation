function intensity = water_drop_scattering(a, lambda, theta, varargin)
p = inputParser;
p.addRequired('a', @(x) isnumeric(x) && isscalar(x));
p.addRequired('lambda', @(x) isnumeric(x) && isvector(x));
p.addRequired('theta', @(x) isnumeric(x) && isvector(x));
p.addParameter('SunSize', 0.5, @(x) isnumeric(x) && isscalar(x));
p.parse(a, lambda, theta, varargin{:});

m = water_refractive_index(lambda);

dq = abs(lambda(2) - lambda(1));
if p.Results.SunSize > 0
    r = p.Results.SunSize / 2;
    smoothing_kernel = sqrt(1 - ((-r:dq:r) / r).^2);
    smoothing_kernel = smoothing_kernel / sum(smoothing_kernel);
    pad_size = ceil(length(smoothing_kernel) / 2);
end

intensity = zeros(length(theta), length(lambda));
for i = 1:length(lambda)
    fprintf('Computing lambda %.3fum, #%d/%d\n', lambda(i), i, length(lambda));
    curr_intensity = mie_theory_scattering(a, m(i), lambda(i), theta);
    curr_intensity = curr_intensity / sum(curr_intensity);
    if p.Results.SunSize > 0
        curr_intensity = padarray(curr_intensity, [1, 0] * pad_size, 'both', 'replicate');
        curr_intensity = conv(curr_intensity, smoothing_kernel, 'same');
        curr_intensity = curr_intensity(pad_size+1:end-pad_size);
    end
    intensity(:, i) = curr_intensity;
end
end