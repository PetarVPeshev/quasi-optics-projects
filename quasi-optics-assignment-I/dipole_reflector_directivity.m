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
h_sweep.Hmin = 0e-3;
h_sweep.Hmax = 25e-3;
h_sweep.Nh = 500;
Nphi = 1200;
Ntheta = 300;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
antenna.W = wave.wavelength / antenna.W_ratio;
antenna.L = wave.wavelength / antenna.L_ratio;
antenna.h = linspace(h_sweep.Hmin, h_sweep.Hmax, h_sweep.Nh);

%% SPHERICAL COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, ~] = wave_vector(medium.er, wave.k0, sph_grid);

%% CALCULATE DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k0, k_comp, 'E', 'J');

%% PERFORM HEIGHT SWEEP
dir_vector = NaN(1, length(antenna.h));
for idx = 1 : 1 : length(antenna.h)

    %% CALCULATE CURRENT FOURIER TRANSFORM
    Jft = ft_current(wave.k0, k_comp, antenna.W, antenna.L, ...
        'dipole', 'x') .* 2j .* sin(k_comp(:, :, 3) * antenna.h(idx));
    
    %% CALCULATE ELECTRIC FAR-FIELD
    Efar = farfield(wave.k0, R, sph_grid, k_comp(:, :, 3), SGFej, Jft);

    %% DIRECTIVITY
    [dir, ~, ~] = directivity(medium.er, Efar, sph_grid, R);
    dir_vector(idx) = dir(1, 1);

end

%% PLOT DIRECTIVITY AS FUNCTION OF HEIGHT IN DECIBEL
figure('Position', [250 250 700 400]);
plot(antenna.h * 1e3, 10*log10(dir_vector), 'LineWidth', 3.0);
grid on;
xlim([0 25]);
xlabel('h / mm');
ylabel('D(\theta=0^{\circ},\phi=0^{\circ}) / dB');
title('Dipole Broadside Directivity');

%% PLOT DIRECTIVITY AS FUNCTION OF HEIGHT IN LINEAR
figure('Position', [250 250 700 400]);
plot(antenna.h * 1e3, dir_vector, 'LineWidth', 3.0);
grid on;
xlim([0 25]);
xlabel('h / mm');
ylabel('D(\theta=0^{\circ},\phi=0^{\circ}) / dB');
title('Dipole Broadside Directivity');
