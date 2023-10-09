close all;
clear;
clc;

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
medium.er = 1;
wave.f = 28e9;
antenna.W_ratio = 40;
antenna.L_ratio = 2;
antenna.H = 5.36e-3;
Nphi = 1200;
Ntheta = 300;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
antenna.W = wave.wavelength / antenna.W_ratio;
antenna.L = wave.wavelength / antenna.L_ratio;

%% SPHERICAL COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, ~] = wave_vector(medium.er, wave.k0, sph_grid);

%% CALCULATE DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k0, k_comp, 'E', 'J');

%% CALCULATE CURRENT FOURIER TRANSFORM
Jft = ft_current(wave.k0, k_comp, antenna.W, antenna.L, 'dipole', 'x') ...
    .* 2j .* sin(k_comp(:, :, 3) * antenna.H);

%% CALCULATE ELECTRIC FAR-FIELD
Efar = farfield(wave.k0, R, sph_grid, k_comp(:, :, 3), SGFej, Jft);
Efar_total = total_field(Efar);

%% UV COORDINATES
uv_mesh = uv_repr(sph_grid);

%% PLOT ELECTRIC FAR-FIELD USING UV REPRESENTATION VIEW (-37.5, 30)
figure();
surface(uv_mesh(:, :, 1), uv_mesh(:, :, 2), norm_magnitude(Efar_total), ...
    'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
zlim([-10 0]);
caxis([-10 0]);
view(-37.5, 30);
xlabel('U');
ylabel('V');
zlabel('|E_{t}| / dB');
title(['Dipole With Reflector Total E-Farfield @ h = ' ...
    num2str(antenna.H * 1e3) ' mm']);

%% PLOT ELECTRIC FAR-FIELD USING UV REPRESENTATION VIEW (0, 90)
figure();
surface(uv_mesh(:, :, 1), uv_mesh(:, :, 2), norm_magnitude(Efar_total), ...
    'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
zlim([-10 0]);
caxis([-10 0]);
view(0, 90);
xlabel('U');
ylabel('V');
zlabel('|E_{t}| / dB');
title(['Dipole With Reflector Total E-Farfield @ h = ' ...
    num2str(antenna.H * 1e3) ' mm']);

%% PLOT ELECTRIC FAR-FIELD IN PHI 0, 45, 90 DEGREES
figure('Position', [250 250 700 400]);
plot(theta * 180 / pi, norm_magnitude(Efar_total(1, :), 'dB'), ...
    'LineWidth', 3.0, 'DisplayName', '\phi = 0^{\circ}');
hold on;
plot(theta * 180 / pi, norm_magnitude(Efar_total(151, :), 'dB'), ...
    'LineWidth', 3.0, 'DisplayName', '\phi = 45^{\circ}');
hold on;
plot(theta * 180 / pi, norm_magnitude(Efar_total(301, :), 'dB'), ...
    'LineWidth', 3.0, 'DisplayName', '\phi = 90^{\circ}');
grid on;
xlim([0 90]);
ylim([-20 0]);
xticks(0 : 10 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E_{t}| / dB');
title(['Dipole With Reflector Total E-Field @ h = ' ...
    num2str(antenna.H * 1e3) ' mm']);
