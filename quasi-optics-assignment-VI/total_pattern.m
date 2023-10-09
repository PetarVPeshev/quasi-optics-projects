close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
wave.f = 10e9;
medium.er = 1;
fss.Nx = 10;
fss.Ny = 10;
dipole.l = 15e-3;
dipole.w = 1e-3;
fss.dx = 20e-3;
fss.dy = 20e-3;
theta_i = [0 30] * pi / 180;
phi_i = [0 0] * pi / 180;
Ntheta = 1200;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;

%% SPHERICAL COORDINATE GRID
theta = linspace(eps, pi / 2, Ntheta);
phi = [0 90 180 270] * pi / 180;
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, k] = wave_vector(medium.er, wave.k, sph_grid);

Efar = zeros( [size(sph_grid, 1, 2), 3, length(theta_i)] );
Efar_total = NaN( [size(sph_grid, 1, 2), length(theta_i)] );
Efar_norm = NaN( [size(sph_grid, 1, 2), length(theta_i)] );
for idx = 1 : 1 : length(theta_i)
    %% INCIDENT SPHERICAL GRID
    sph_grid_i = NaN(1, 1, 2);
    sph_grid_i(:, :, 1) = theta_i(idx);
    sph_grid_i(:, :, 2) = phi_i(idx);

    %% ARRAY FACTOR
    AF = array_factor(fss.Nx, fss.Ny, fss.dx, fss.dy, k, k_comp, ...
        theta_i(idx), phi_i(idx));

    %% BF CURRENT
    ibf = bf_current(medium.er, dipole.l, dipole.w, fss.dx, fss.dy, ...
        wave.k, sph_grid_i);

    %% WINDOWED CURRENT SPECTRUM
    Jw = windowed_current(dipole.l, dipole.w, ibf, AF, k, k_comp);

    %% DYADIC SPECTRAL GREEN'S FUNCTION
    SGFej = dyadic_sgf(medium.er, wave.k, k_comp, 'E', 'J');

    %% ACTIVE ELEMENT ELECTRIC FAR-FIELD
    Efar(:, :, 1, idx) = 1j * k_comp(:, :, 3) .* SGFej(:, :, 1, 1) ...
        .* Jw(:, :, 1) .* exp(-1j * k * R) ./ (2 * pi * R);
    Efar(:, :, 2, idx) = 1j * k_comp(:, :, 3) .* SGFej(:, :, 2, 1) ...
        .* Jw(:, :, 1) .* exp(-1j * k * R) ./ (2 * pi * R);
    Efar(:, :, 3, idx) = 1j * k_comp(:, :, 3) .* SGFej(:, :, 3, 1) ...
        .* Jw(:, :, 1) .* exp(-1j * k * R) ./ (2 * pi * R);
    Efar_total(:, :, idx) = total_field(Efar(:, :, :, idx));
    Efar_norm(:, :, idx) = norm_magnitude(Efar_total(:, :, idx), 'dB');
end

%% PLOT E AND H PLANE IMPEDANCE
theta_length = length(theta);
theta_plot = NaN(1, theta_length * 2);
theta_plot(1 : theta_length) = - fliplr(theta) * 180 / pi;
theta_plot(theta_length + 1 : end) = theta * 180 / pi;
figure('Position', [250 250 900 400]);
for idx = 1 : 1 : length(theta_i)
    plane_field = NaN(1, theta_length * 2);
    plane_field(1 : theta_length) = fliplr(Efar_norm(3, :, idx));
    plane_field(theta_length + 1 : end) = Efar_norm(1, :, idx);
    plot(theta_plot, plane_field, 'LineWidth', 2.0, ...
        'DisplayName', ['|E|, \theta_{i} = ' ...
        num2str(theta_i(idx) * 180 / pi) ' deg and \phi_{i} = ' ...
        num2str(phi_i(idx) * 180 / pi) ' deg']);
    hold on;
end
hold off;
grid on;
legend show;
legend('location', 'bestoutside');
xticks(-90 : 15 : 90);
xlim([-90 90]);
ylim([-30 0]);
xlabel('\theta / deg');
ylabel('|E| / dB');
title(['E Far-Field @ E-plane, f = ' num2str(wave.f * 1e-9) ...
    ' GHz, L = ' num2str(dipole.l * 1e3) ' mm, and W = ' ...
    num2str(dipole.w * 1e3) ' mm, d_{x} = ' num2str(fss.dx * 1e3) ...
    ' mm, and d_{y} = ' num2str(fss.dy * 1e3) ' mm']);
saveas(gcf, 'figures\E_plane.fig');
