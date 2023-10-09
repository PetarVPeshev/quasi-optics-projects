close all;
clear;
clc;

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
medium.er = 1;
fsweep.fmin = 15e9;
fsweep.fmax = 195e9;
fsweep.Nf = 180;
antenna.W = 0.268e-3;
antenna.L = 5.36e-3;
Nphi = 1200;
Ntheta = 300;
R = 1;

%% DEPENDENT PARAMETERS
wave.f = linspace(fsweep.fmin, fsweep.fmax, fsweep.Nf);
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;

%% SPHERICAL COORDINATE GRID
theta = linspace(eps, pi, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% PERFORM FREQUENCY SWEEP
dir_vector = NaN(1, length(wave.f));
for idx = 1 : 1 : length(wave.f)

    %% WAVE VECTOR COMPONENTS
    [k_comp, ~] = wave_vector(medium.er, wave.k0(idx), sph_grid);

    %% CALCULATE DYADIC SPECTRAL GREEN'S FUNCTION
    SGFej = dyadic_sgf(medium.er, wave.k0(idx), k_comp, 'E', 'J');
    
    %% CALCULATE CURRENT FOURIER TRANSFORM
    Jft = ft_current(wave.k0(idx), k_comp, antenna.W, antenna.L, ...
        'dipole', 'x');
    
    %% CALCULATE ELECTRIC FAR-FIELD
    Efar = farfield(wave.k0(idx), R, sph_grid, k_comp(:, :, 3), SGFej, ...
        Jft);

    %% DIRECTIVITY
    [dir, ~, ~] = directivity(medium.er, Efar, sph_grid, R);
    dir_vector(idx) = dir(1, 1);

end

%% PLOT DIRECTIVITY AS FUNCTION OF FREQUENCY
figure('Position', [250 250 700 400]);
plot(wave.f * 1e-9, 10*log10(dir_vector), 'LineWidth', 3.0);
grid on;
xlim([15 195]);
xlabel('f / GHz');
ylabel('D(\theta=0^{\circ},\phi=0^{\circ}) / dB');
title('Dipole Broadside Directivity');
