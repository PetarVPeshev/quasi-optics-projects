close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
wave.f = 50e9;
medium.er = 1;
reflector.f = 3;
Nphi = 2000;
Ntheta = 2000;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;
feed.D = 4 * wave.wavelength;
reflector.D = reflector.f / 0.6;

%% COORDINATES GRID
theta = linspace(eps, pi / 2 - eps, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, ~] = wave_vector(medium.er, wave.k, sph_grid);

%% DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k, k_comp, 'E', 'J');

%% CIRCULAR FEED CURRENT DENSITY FOURIER TRANSFORM
feed.Jft = ft_current(wave.k, feed.D / 2, sph_grid(:, :, 1), ...
    'circular', 'y');

%% ELECTRIC FAR-FIELD OF CIRCULAR FEED
feed.Efar_comp = farfield(wave.k, R, sph_grid, k_comp(:, :, 3), ...
    SGFej, feed.Jft) * 2 * pi * R / exp(-1j * wave.k * R);

%% APERTURE CURRENT
[~, reflector.M, cart_grid] = aperture_current(feed.Efar_comp, wave.k, ...
    medium.er, sph_grid, 'reflector', reflector.f, reflector.D, 'Love');

%% PLOT MAGNETIC CURRENT DENSITY
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(cart_grid(:, :, 1), cart_grid(:, :, 2), ...
    norm_magnitude(reflector.M(:, :, 1), 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([-reflector.D reflector.D] / 2);
ylim([-reflector.D reflector.D] / 2);
caxis([-40 0]);
view(0, 90);
xlabel('x / m');
ylabel('y / m');
zlabel('|M_{x}| / dB');
title('M_{x}');
subplot(1, 2, 2);
surface(cart_grid(:, :, 1), cart_grid(:, :, 2), ...
    norm_magnitude(reflector.M(:, :, 2), 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([-reflector.D reflector.D] / 2);
ylim([-reflector.D reflector.D] / 2);
caxis([-40 0]);
view(0, 90);
xlabel('x / m');
ylabel('y / m');
zlabel('|M_{y}| / dB');
title('M_{y}');
sgtitle(['Aperture Magnetic Current Density @  f = ' ...
    num2str(wave.f * 1e-9) ' GHz, D_{feed}=4\lambda, F = ' ...
    num2str(reflector.f) ' m, and D = ' num2str(reflector.D) ' m'], ...
    'FontSize', 17, 'FontWeight', 'bold');
saveas(gcf, 'figures\aperture_M.fig');
