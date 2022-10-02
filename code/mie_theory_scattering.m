function [intensity, Q_sct, Q_ext] = mie_theory_scattering(drop_radius, refractive_index, wavelength, theta, varargin)
% Ricatti Bessel function of first and second kind
%   psi_n(z) = sqrt(pi * z / 2) * J_{n+1/2}(z)
%   chi_n(z) = (-1)^n * sqrt(pi * z / 2) * J_{-(n+1/2)}(z)
%
% They are related to SphericalBesselJ (spherical Bessel function of the first kind) and
% SphericalBesselY (spherical Bessel function of the second kind) function in mathematica
%   j_n(z) = sqrt(pi / 2 / z) * J_{n+1/2}(z)
%   y_n(z) = sqrt(pi / 2 / z) * Y_{n+1/2}(z)
% thus,
%   psi_n(z) = z * j_n(z)
%   chi_n(z) = -z * y_n(z)
%
% Ricatti Hankel function:
%   zeta_n(z) = psi_n(z) + i * chi_n(z)
%
% It is related to SphericalHankelH2 (spherical Hankel function of the second kind) function
% in mathematica
%   h2_n(z) = j_n(z) - i * y_n(z) = sqrt(pi / 2 / z) * H2_{n+1/2}(z)
% thus,
%   zeta_n(z) = z * j_n(z) - i * z * y_n(z) = z * h2_n(z)

x = 2 * pi ./ wavelength * drop_radius;
y = refractive_index * x;

if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1}, 'direct')
    [intensity, Q_sct, Q_ext] = direct_method(x, y, refractive_index, theta);
else
    [intensity, Q_sct, Q_ext] = iterate_method(x, y, refractive_index, theta);
end
end


function [intensity, Q_sct, Q_ext] = direct_method(x, y, m, theta)
psi = @(n, x) sqrt(pi * x / 2) .* besselj(n + 0.5, x);
zeta = @(n, x) sqrt(pi * x / 2) .* besselh(n + 0.5, 2, x);

dpsi = @(n, x) psi(n-1, x) - psi(n, x) * n ./ x;
dzeta = @(n, x) zeta(n-1, x) - zeta(n, x) * n ./ x;

theta = theta(:) * pi / 180;        % Convert to radian
mu = cos(theta);

P_n_prev = ones(size(mu));
P_n_curr = mu;

S1 = 0;
S2 = 0;
ab0 = [];
Q_sct = 0;
Q_ext = 0;
N_MAX = 300000;
for n = 1:N_MAX
    pi_n = n ./ (mu.^2 - 1) .* (mu .* P_n_curr - P_n_prev);
    tau_n = 1 ./ (mu.^2 - 1) .* (n * (1 + n - n * mu) .* P_n_curr - (n * mu) .* P_n_prev);

    an = (dpsi(n, y) .* psi(n, x) - m * psi(n, y) .* dpsi(n, x)) ./ ...
        (dpsi(n, y) .* zeta(n, x) - m * psi(n, y) .* dzeta(n, x));
    bn = (m * dpsi(n, y) .* psi(n, x) - psi(n, y) .* dpsi(n, x)) ./ ...
        (m * dpsi(n, y) .* zeta(n, x) - psi(n, y) .* dzeta(n, x));
    ab = abs(an).^2 + abs(bn).^2;
    Q_sct = Q_sct + (2 * n + 1) * ab;
    Q_ext = Q_ext + (2 * n + 1) * real(an + bn);

    if isempty(ab0)
        ab0 = ab;
    end
    if all(ab ./ ab0 < 1e-14)
        break;
    end

    S1 = S1 + (2 * n + 1) / (n * (n + 1)) .* (an .* pi_n + bn .* tau_n);
    S2 = S2 + (2 * n + 1) / (n * (n + 1)) .* (bn .* pi_n + an .* tau_n);

    p = legendre(n + 1, mu);
    P_n_prev = P_n_curr;
    P_n_curr = p(1, :)';
end
intensity = (abs(S1).^2 + abs(S2).^2);
Q_sct = Q_sct * 2 ./ x.^2;
Q_ext = Q_ext * 2 ./ x.^2;
end


function [intensity, Q_sct, Q_ext] = iterate_method(x, y, m, theta)
theta = theta(:) * pi / 180;        % Convert to radian
mu = cos(theta);

psi_x_prev = sqrt(pi * x / 2) .* besselj(0.5, x);
psi_x_curr = sqrt(pi * x / 2) .* besselj(1.5, x);
zeta_x_prev = sqrt(pi * x / 2) .* besselh(0.5, 2, x);
zeta_x_curr = sqrt(pi * x / 2) .* besselh(1.5, 2, x);
psi_y_prev = sqrt(pi * y / 2) .* besselj(0.5, y);
psi_y_curr = sqrt(pi * y / 2) .* besselj(1.5, y);

pi_n_prev = zeros(size(mu));
pi_n_curr = ones(size(mu));

S1_real = 0;
S1_imag = 0;
S2_real = 0;
S2_imag = 0;
ab0 = [];
Q_sct = 0;
Q_ext = 0;
N_MAX = 300000;
% N_MAX = ceil(x + 4 * nthroot(x, 3) + 2);
for n = 1:N_MAX
    tau_n_curr = n * mu .* pi_n_curr - (n + 1) * pi_n_prev;

    dpsi_x = -n ./ x .* psi_x_curr + psi_x_prev;
    dzeta_x = -n ./ x .* zeta_x_curr + zeta_x_prev;
    dpsi_y = -n ./ y .* psi_y_curr + psi_y_prev;

    an = (dpsi_y .* psi_x_curr - m * psi_y_curr .* dpsi_x) ./ ...
        (dpsi_y .* zeta_x_curr - m * psi_y_curr .* dzeta_x);
    bn = (m * dpsi_y .* psi_x_curr - psi_y_curr .* dpsi_x) ./ ...
        (m * dpsi_y .* zeta_x_curr - psi_y_curr .* dzeta_x);
    ab = abs(an).^2 + abs(bn).^2;
    Q_sct = Q_sct + (2 * n + 1) * ab;
    Q_ext = Q_ext + (2 * n + 1) * real(an + bn);

    if isempty(ab0)
        ab0 = ab;
    end
    if all(ab ./ ab0 < 1e-14)
        break;
    end

    k = (2 * n + 1) / (n * (n + 1));
    an_real = real(an);
    an_imag = imag(an);
    bn_real = real(bn);
    bn_imag = imag(bn);
    S1_real = S1_real + k .* (an_real .* pi_n_curr + bn_real .* tau_n_curr);
    S1_imag = S1_imag + k .* (an_imag .* pi_n_curr + bn_imag .* tau_n_curr);
    S2_real = S2_real + k .* (bn_real .* pi_n_curr + an_real .* tau_n_curr);
    S2_imag = S2_imag + k .* (bn_imag .* pi_n_curr + an_imag .* tau_n_curr);

    pi_n = (2 * n + 1) / n * mu .* pi_n_curr - (n + 1) / n * pi_n_prev;
    pi_n_prev = pi_n_curr;
    pi_n_curr = pi_n;

    psi_x = (2 * n + 1) ./ x .* psi_x_curr - psi_x_prev;
    psi_x_prev = psi_x_curr;
    psi_x_curr = psi_x;

    psi_y = (2 * n + 1) ./ y .* psi_y_curr - psi_y_prev;
    psi_y_prev = psi_y_curr;
    psi_y_curr = psi_y;

    zeta_x = (2 * n + 1) ./ x .* zeta_x_curr - zeta_x_prev;
    zeta_x_prev = zeta_x_curr;
    zeta_x_curr = zeta_x;
end
S1 = S1_real + 1i * S1_imag;
S2 = S2_real + 1i * S2_imag;
intensity = (abs(S1).^2 + abs(S2).^2);
Q_sct = Q_sct * 2 ./ x.^2;
Q_ext = Q_ext * 2 ./ x.^2;
end