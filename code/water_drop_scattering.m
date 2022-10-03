function intensity = water_drop_scattering(a, lambda, theta, varargin)
p = inputParser;
p.addRequired('a', @(x) isnumeric(x) && isscalar(x));
p.addRequired('lambda', @(x) isnumeric(x) && isvector(x));
p.addRequired('theta', @(x) isnumeric(x) && isvector(x));
p.addParameter('SunSize', 0.5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Parallel', false, @(x) islogical(x) && isscalar(x));
p.parse(a, lambda, theta, varargin{:});

m = water_refractive_index(lambda);

dq = abs(theta(2) - theta(1));
if p.Results.SunSize > 0
    r = p.Results.SunSize / 2;
    smoothing_kernel = sqrt(1 - ((-r:dq:r) / r).^2);
    smoothing_kernel = smoothing_kernel / sum(smoothing_kernel);
end

intensity = zeros(length(theta), length(lambda));
if p.Results.Parallel
    parfor i = 1:length(lambda)
        curr_intensity = mie_theory_scattering(a, m(i), lambda(i), theta);
        intensity(:, i) = curr_intensity;
    end
else
    for i = 1:length(lambda)
        curr_intensity = mie_theory_scattering(a, m(i), lambda(i), theta);
        intensity(:, i) = curr_intensity;
    end
end
if p.Results.SunSize > 0
    intensity = imfilter(intensity, smoothing_kernel(:), 'same', 'symmetric');
end
end