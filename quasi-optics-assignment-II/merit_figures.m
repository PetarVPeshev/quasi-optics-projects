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
Nphi = 4000;
Ntheta = 4000;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;
feed.D = 4 * wave.wavelength;
reflector.D = reflector.f ./ ( 0.6 : 0.05 : 6 );

%% COORDINATES GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, ~] = wave_vector(medium.er, wave.k, sph_grid);

%% DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k, k_comp, 'E', 'J');

%% FIGURES OF MERIT
eta_tap = NaN(1, length(reflector.D));
eta_so = NaN(1, length(reflector.D));
dir_max = NaN(1, length(reflector.D));

for idx = 1 : 1 : length(reflector.D)

    %% MAXIMUM DIRECTIVITY
    dir_max(idx) = 4 * ( ( pi * reflector.D(idx) ) ^ 2 ) ...
        / ( 4 * ( wave.wavelength ^ 2 ) );

    %% CIRCULAR FEED CURRENT DENSITY FOURIER TRANSFORM
    Jft_feed = ft_current(wave.k, feed.D / 2, sph_grid(:, :, 1), ...
        'circular', 'y');

    %% ELECTRIC FAR-FIELD OF CIRCULAR FEED
    Efar_feed = farfield(wave.k, R, sph_grid, k_comp(:, :, 3), ...
        SGFej, Jft_feed);
    Efar_comp_feed = Efar_feed * 2 * pi * R / exp(-1j * wave.k * R);

    %% TAPER EFFICIENCY
    [eta_tap(idx), ~] = taper_efficiency(Efar_comp_feed, wave.k, ...
        sph_grid, 'reflector', reflector.f, reflector.D(idx));

    %% SPILLOVER EFFICIENCY
    eta_so(idx) = spillover_efficiency(Efar_feed, medium.er, sph_grid, ...
        R, 'reflector', reflector.f, reflector.D(idx));

end
eta_tap(eta_tap > 1) = 1;

%% APERTURE EFFICIENCY
eta_ap = eta_tap .* eta_so;

%% DIRECTIVITY
dir = dir_max .* eta_tap;

%% GAIN
gain = dir_max .* eta_ap;

%% PLOT TAPER, SPILLOVER, AND APERTURE EFFICIENCY
figure('Position', [250 250 725 400]);
plot(reflector.D, eta_tap, 'LineWidth', 3.0, 'DisplayName', '\eta_{tap}');
hold on;
plot(reflector.D, eta_so, 'LineWidth', 3.0, 'DisplayName', '\eta_{so}');
hold on;
plot(reflector.D, eta_ap, 'LineWidth', 3.0, 'DisplayName', '\eta_{ap}');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('D / m');
ylabel('\eta');
title(['Efficiency @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    'D_{feed} = 4\lambda, and F = ' num2str(reflector.f) ' m']);
saveas(gcf, 'figures\efficiency.fig');

%% PLOT MAXIMUM DIRECTIVITY, DIRECTIVITY, AND GAIN
figure('Position', [250 250 725 400]);
plot(reflector.D, 10 * log10(dir_max), 'LineWidth', 3.0, ...
    'DisplayName', 'max directivity');
hold on;
plot(reflector.D, 10 * log10(dir), 'LineWidth', 3.0, ...
    'DisplayName', 'directivity');
hold on;
plot(reflector.D, 10 * log10(gain), 'LineWidth', 3.0, ...
    'DisplayName', 'gain');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('D / m');
ylabel('- / dB');
title(['Directivity and Gain @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz, D_{feed} = 4\lambda, and F = ' num2str(reflector.f) ' m']);
saveas(gcf, 'figures\directivity_and_gain.fig');
