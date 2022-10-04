function [intensity, lambda] = water_drop_scattering(a, lambda, theta, varargin)
p = inputParser;
p.addRequired('a', @(x) isnumeric(x) && isscalar(x));
p.addRequired('lambda', @(x) isnumeric(x) && isvector(x));
p.addRequired('theta', @(x) isnumeric(x) && isvector(x));
p.addParameter('SunSize', 0.5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('LambdaRelErr', 0.1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('AdaptiveLambda', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Parallel', false, @(x) islogical(x) && isscalar(x));
p.parse(a, lambda, theta, varargin{:});

dq = abs(theta(2) - theta(1));
if p.Results.SunSize > 0
    r = p.Results.SunSize / 2;
    smoothing_kernel = sqrt(1 - ((-r:dq:r) / r).^2);
    smoothing_kernel = smoothing_kernel / sum(smoothing_kernel);
end

if p.Results.AdaptiveLambda
    [intensity, lambda] = adaptive_sample_lambda(a, lambda, theta, p.Results.LambdaRelErr, p.Results.Parallel);
else
    intensity = zeros(length(theta), length(lambda));
    m = water_refractive_index(lambda);
    if p.Results.Parallel
        parfor i = 1:length(lambda)
            intensity(:, i) = mie_scattering(a, m(i), lambda(i), theta);
        end
    else
        for i = 1:length(lambda)
            intensity(:, i) = mie_scattering(a, m(i), lambda(i), theta);
        end
    end
end
if p.Results.SunSize > 0
    intensity = imfilter(intensity, smoothing_kernel(:), 'same', 'symmetric');
end
end


function [intensity, lambda] = adaptive_sample_lambda(a, lambda_lim, theta, rel_err_lim, parallel)
lambda0 = linspace(lambda_lim(1), lambda_lim(2), 50);
intensity0 = zeros(length(theta), length(lambda0));
if ~parallel
    for i = 1:length(lambda0)
        m = water_refractive_index(lambda0(i));
        intensity0(:, i) = mie_scattering(a, m, lambda0(i), theta);
    end
else
    parfor i = 1:length(lambda0)
        m = water_refractive_index(lambda0(i));
        intensity0(:, i) = mie_scattering(a, m, lambda0(i), theta);
    end
end
    
rel_err = inf;
dw = lambda0(2) - lambda0(1);
while rel_err > rel_err_lim && dw > 0.0001
    lambda1 = (lambda0(1:end-1) + lambda0(2:end)) / 2;
    intensity1 = zeros(length(theta), length(lambda1));
    if ~parallel
        for i = 1:length(lambda1)
            m = water_refractive_index(lambda1(i));
            intensity1(:, i) = mie_scattering(a, m, lambda1(i), theta);
        end
    else
        parfor i = 1:length(lambda1)
            m = water_refractive_index(lambda1(i));
            intensity1(:, i) = mie_scattering(a, m, lambda1(i), theta);
        end
    end

    err = (intensity0(:, 1:end-1) + intensity0(:, 2:end)) / 2 - intensity1;
    rel_err = sum(abs(err(:))) ./ sum(abs(intensity1(:)));

    lambda = zeros(1, length(lambda0) + length(lambda1));
    lambda(1:2:end) = lambda0;
    lambda(2:2:end) = lambda1;
    lambda0 = lambda;
    dw = lambda0(2) - lambda0(1);
    
    intensity = zeros(length(theta), length(lambda0));
    intensity(:, 1:2:end) = intensity0;
    intensity(:, 2:2:end) = intensity1;
    intensity0 = intensity;
end
end