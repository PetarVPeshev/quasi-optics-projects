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
phi = [0; 45; 90] * pi / 180;
Ntheta = 300;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
antenna.W = wave.wavelength / antenna.W_ratio;
antenna.L = wave.wavelength / antenna.L_ratio;

%% SPHERICAL COORDINATE GRID
theta = linspace(eps, pi, Ntheta);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, ~] = wave_vector(medium.er, wave.k0, sph_grid);

%% CALCULATE DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k0, k_comp, 'E', 'J');

%% CALCULATE CURRENT FOURIER TRANSFORM
Jft = ft_current(wave.k0, k_comp, antenna.W, antenna.L, 'dipole', 'x');

%% CALCULATE ELECTRIC FAR-FIELD
Efar = farfield(wave.k0, R, sph_grid, k_comp(:, :, 3), SGFej, Jft);
Efar_total = total_field(Efar);

%% PLOT TOTAL ELECTRIC FAR-FIELD
figure('Position', [250 250 700 400]);
for idx = 1 : 1 : length(phi)
    plot(theta * 180 / pi, norm_magnitude(Efar_total(idx, :), 'dB'), ...
        'LineWidth', 3.0, 'DisplayName', ['\phi = ' ...
        num2str(phi(idx) * 180 / pi) '^{\circ}']);
    hold on;
end
hold off;
grid on;
xlim([0 180]);
ylim([-20 0]);
xticks(0 : 10 : 180);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E_{t}| / dB');
title('Dipole Total E-Field');
