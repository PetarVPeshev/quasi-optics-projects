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

%% ARRAY FACTOR
% For AEP, AFx = 1 and AFy = 1, envelop of all patterns in the array
AF = ones( [size(sph_grid, 1, 2), 2] );

%% BF CURRENT
ibf = bf_current(medium.er, dipole.l, dipole.w, fss.dx, fss.dy, ...
    wave.k, sph_grid);

%% WINDOWED CURRENT SPECTRUM
Jw = windowed_current(dipole.l, dipole.w, ibf, AF, k, k_comp);

%% DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k, k_comp, 'E', 'J');

%% ACTIVE ELEMENT ELECTRIC FAR-FIELD
Efar = zeros( [size(sph_grid, 1, 2), 3] );
Efar(:, :, 1) = 1j * k_comp(:, :, 3) .* SGFej(:, :, 1, 1) ...
    .* Jw(:, :, 1) .* exp(-1j * k * R) ./ (2 * pi * R);
Efar(:, :, 2) = 1j * k_comp(:, :, 3) .* SGFej(:, :, 2, 1) ...
    .* Jw(:, :, 1) .* exp(-1j * k * R) ./ (2 * pi * R);
Efar(:, :, 3) = 1j * k_comp(:, :, 3) .* SGFej(:, :, 3, 1) ...
    .* Jw(:, :, 1) .* exp(-1j * k * R) ./ (2 * pi * R);
Efar_total = total_field(Efar);

%% FLOQUET MODES
[~, ~, modes] = floquet_modes(wave.k, medium.er, fss.dx, fss.dy, sph_grid);

%% GRATING LOBES
lobes = grating_lobes(k, fss.dx, fss.dy, modes);

%% PLOT E AND H PLANE IMPEDANCE
Efar_norm = norm_magnitude(Efar_total, 'dB');
theta_length = length(theta);
theta_plot = NaN(1, theta_length * 2);
theta_plot(1 : theta_length) = - fliplr(theta) * 180 / pi;
theta_plot(theta_length + 1 : end) = theta * 180 / pi;
labels = ["E-plane", "H-plane"];
figure('Position', [250 250 800 400]);
for idx = 1 : 1 : 2
    plane_field = NaN(1, theta_length * 2);
    plane_field(1 : theta_length) = fliplr(Efar_norm(idx + 2, :));
    plane_field(theta_length + 1 : end) = Efar_norm(idx, :);
    plot(theta_plot, plane_field, 'LineWidth', 2.0, ...
        'DisplayName', [labels(idx)]);
    hold on;
end
hold off;
grid on;
legend show;
legend('location', 'bestoutside');
xticks(-90 : 15 : 90);
xlim([-90 90]);
ylim([-25 0]);
xlabel('\theta / deg');
ylabel('E / dB');
title(['Active Element Pattern @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz, L = ' num2str(dipole.l * 1e3) ' mm, W = ' ...
    num2str(dipole.w * 1e3) ' mm, d_{x} = ' num2str(fss.dx * 1e3) ...
    ' mm, and d_{y} = ' num2str(fss.dy * 1e3) ' mm']);
saveas(gcf, 'figures\aep.fig');

%% PLOT GRATING LOBES
planes_theta = linspace(- pi / 2, pi / 2, 2001);
planes_k = NaN( [2, size(planes_theta)] );
planes_k(1, :) = wave.k * sin(planes_theta);
planes_k(2, :) = 0;
figure('Position', [250 250 600 400]);
plot(lobes(:, :, 1)', lobes(:, :, 2)', ...
    'Color', [0 0 0], 'LineWidth', 2.0);
hold on;
patch(lobes(5, :, 1), lobes(5, :, 2), 'black', ...
    'FaceAlpha', 0.1);
hold on;
plot(planes_k(1, :), planes_k(2, :), ...
    'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2.0);
hold on;
plot(planes_k(2, :), planes_k(1, :), ...
    'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0);
grid on;
legend show;
legend('lobes', '', '', '', '', '', '', '', '', 'visible region', ...
    'E-plane', 'H-plane', 'location', 'bestoutside');
set(gca, 'XTickLabel', [], 'YTickLabel', []);
xlabel('k_{mx}');
ylabel('k_{my}');
title(['Grating Lobes @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    'd_{x} = ' num2str(fss.dx * 1e3) ' mm, and d_{y} = ' ...
    num2str(fss.dy * 1e3) ' mm']);
saveas(gcf, 'figures\grating_lobes.fig');
