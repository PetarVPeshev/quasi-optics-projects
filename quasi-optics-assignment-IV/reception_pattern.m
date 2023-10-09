close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

if ~exist([pwd() '\results'], 'dir')
    mkdir('results');
end

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
wave.f = 50e9;
wave.E0 = 1;
medium.er = 1;
theta_ib = [0 3] * pi / 180;
phi_ib = [0 360] * pi / 180;
Nphi = 200;
Ntheta = 200;
Nphi_i = 200;
Ntheta_i = 100;
R = 1;
planes = [0, 90];

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;
feed.D = 4 * wave.wavelength;
reflector.D = 100 * wave.wavelength;
reflector.f = reflector.D * 2;

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% INCIDENT ANGLE COORDINATE GRID
theta_i = linspace(theta_ib(1), theta_ib(2), Ntheta_i);
phi_i = linspace(phi_ib(1), phi_ib(2), Nphi_i);
sph_grid_i = meshgrid_comb(theta_i, phi_i);

%% TRANSMISSION A FIELD
[tx.Va, tx.rad_power] = feed_va_field(feed.D / 2, wave.k, sph_grid);
tx.Va_total = total_field(tx.Va);

%% RECEPTION GO FIELD
[rx.Ego, rx.Vgo, rx.sph_grid] = go_field(reflector.f, reflector.D, ...
    wave.E0, wave.k, sph_grid, sph_grid_i);
rx.Vgo_total = total_field(rx.Vgo, '5D');

%% RECEPTION POWER AND EFFICIENCY
[rx.Prx, rx.eta] = rx_power(rx.Vgo, tx.Va, tx.rad_power, wave.E0, ...
    reflector.D / 2, sph_grid, rx.sph_grid);
rx.Prx_norm = norm_magnitude(rx.Prx, 'dB');

%% RECEPTION DIRECTIVITY AND GAIN
rx.dir_max = 4 * ((pi * reflector.D) ^ 2) / (4 * (wave.wavelength ^ 2));
rx.dir = 4 * pi / (sum( sum((rx.Prx / max(rx.Prx, [], 'all')) ...
    .* sin(sph_grid_i(:, :, 1))) ) * (theta_i(2) - theta_i(1)) ...
    * (phi_i(2) - phi_i(1)));
rx.dir_db = 10 * log10(rx.dir);
rx.gain = rx.dir_max * rx.eta;
rx.gain_db = 10 * log10(rx.gain);
rx.gain_max = max(rx.gain, [], 'all');
rx.gain_max_db = 10 * log10(rx.gain_max);

%% TRANSMISSION WAVE VECTOR COMPONENTS
[tx.feed.k_comp, ~] = wave_vector(medium.er, wave.k, sph_grid);
[tx.reflector.k_comp, ~] = wave_vector(medium.er, wave.k, sph_grid_i);

%% TRANSMISSION DYADIC SPECTRAL GREEN'S FUNCTION
tx.feed.SGFej = dyadic_sgf(medium.er, wave.k, tx.feed.k_comp, 'E', 'J');
tx.reflector.SGFej = dyadic_sgf(medium.er, wave.k, tx.reflector.k_comp, ...
    'E', 'J');

%% TRANSMISSION FEED CURRENT DENSITY FOURIER TRANSFORM
tx.feed.Jft = ft_current(wave.k, feed.D / 2, sph_grid(:, :, 1), ...
    'circular', 'y');

%% TRANSMISSION FEED ELECTRIC FAR-FIELD
tx.feed.Ecomp = farfield(wave.k, R, sph_grid, tx.feed.k_comp(:, :, 3), ...
    tx.feed.SGFej, tx.feed.Jft) * 2 * pi * R / exp(-1j * wave.k * R);

%% TRANSMISSION APERTURE CURRENT FOURIER TRANSFORM
[~, ~, tx.reflector.Jft] = aperture_current(tx.feed.Ecomp, wave.k, ...
    medium.er, sph_grid, 'reflector', reflector.f, reflector.D, ...
    'Schelkunoff', 'FT', tx.reflector.k_comp);

%% TRANSMISSION APERTURE ELECTRIC FAR-FIELD
tx.reflector.Efar = farfield(wave.k, R, sph_grid_i, ...
    tx.reflector.k_comp(:, :, 3), tx.reflector.SGFej, tx.reflector.Jft);
tx.reflector.Efar_total = total_field(tx.reflector.Efar);
tx.reflector.Efar_norm = norm_magnitude(tx.reflector.Efar_total, 'dB');

%% TRANSMISSION DIRECTIVITY, EFFICIENCY, AND GAIN
tx.dir_max = rx.dir_max;
[tx.eta_tap, ~] = taper_efficiency(tx.feed.Ecomp, wave.k, ...
    sph_grid, 'reflector', reflector.f, reflector.D);
tx.eta_so = spillover_efficiency(tx.feed.Ecomp, medium.er, sph_grid, ...
    R, 'reflector', reflector.f, reflector.D);
tx.eta_ap = tx.eta_tap * tx.eta_so;
tx.dir = tx.dir_max .* tx.eta_tap;
tx.dir_db = 10 * log10(tx.dir);
tx.gain = tx.dir_max .* tx.eta_ap;
tx.gain_db = 10 * log10(tx.gain);

%% INCIDENT UV COORDINATES
uv_grid_i = uv_repr(sph_grid_i);

%% PLOT RECEPTION PATTERN IN UV REPRESENTATION
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(uv_grid_i(:, :, 1), uv_grid_i(:, :, 2), rx.Prx_norm, ...
    'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
caxis([-40 0]);
view(0, 90);
xlabel('U');
ylabel('V');
zlabel('|P_{rx}| / dBs');
subplot(1, 2, 2);
surface(uv_grid_i(:, :, 1), uv_grid_i(:, :, 2), rx.Prx_norm, ...
    'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
caxis([-40 0]);
view(-37.5, 30);
xlabel('U');
ylabel('V');
zlabel('|P_{rx}| / dB');
sgtitle(['Reception Pattern @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz, D = 100\lambda, D_{feed} = 4\lambda, and f_{#} = 2'], ...
    'FontSize', 17, 'FontWeight', 'bold');
saveas(gcf, 'figures\rx_pattern.fig');

%% PLOT E AND H PLANES RECEPTION AND TRANSMISSION PATTERNS
theta_plot = NaN(1, length(theta_i) * 2);
theta_plot(1 : length(theta_i)) = - fliplr(theta_i) * 180 / pi;
theta_plot(length(theta_i) + 1 : end) = theta_i * 180 / pi;
line_styles = ["-", "--"];
figure('Position', [250 250 725 400]);
for idx = 1 : 1 : length(planes)
    plane_idx = find(round(phi * 180 / pi, 0) == planes(idx), 1);
    plane_field = NaN(1, length(theta_i) * 2);
    plane_field(1 : length(theta_i)) = fliplr(rx.Prx_norm(plane_idx, :));
    plane_field(length(theta_i) + 1 : end) = rx.Prx_norm(plane_idx, :);
    plot(theta_plot, plane_field, [line_styles(idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['P_{rx}, \phi = ' num2str(planes(idx)) ...
        ' deg and ' num2str(planes(idx) + 180) ' deg']);
    hold on;
end
for idx = 1 : 1 : length(planes)
    plane_idx = find(round(phi * 180 / pi, 0) == planes(idx), 1);
    plane_field = NaN(1, length(theta_i) * 2);
    plane_field(1 : length(theta_i)) = ...
        fliplr(tx.reflector.Efar_norm(plane_idx, :));
    plane_field(length(theta_i) + 1 : end) = ...
        tx.reflector.Efar_norm(plane_idx, :);
    plot(theta_plot, plane_field, [line_styles(idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['|E_{tx}|, \phi = ' num2str(planes(idx)) ...
        ' deg and ' num2str(planes(idx) + 180) ' deg']);
    hold on;
end
hold off;
grid on;
ylim([-60 0]);
xlim([-3 3]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|P_{rx}| / dBs and |E_{tx}| / dB');
title(['RX and TX Pattern @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz, D = 100\lambda, D_{feed} = 4\lambda, and f_{#} = 2']);
saveas(gcf, 'figures\rx_and_tx_E_and_H_plane.fig');

%% SAVE RX AND TX STRUCTS
save_var.dir_max = rx.dir_max;
save_var.rx.eta = rx.eta;
save_var.rx.dir = rx.dir;
save_var.rx.dir_db = rx.dir_db;
save_var.rx.gain = rx.gain;
save_var.rx.gain_db = rx.gain_db;
save_var.rx.gain_max = rx.gain_max;
save_var.rx.gain_max_db = rx.gain_max_db;
save_var.tx.eta_tap = tx.eta_tap;
save_var.tx.eta_so = tx.eta_so;
save_var.tx.eta_ap = tx.eta_ap;
save_var.tx.dir = tx.dir;
save_var.tx.dir_db = tx.dir_db;
save_var.tx.gain = tx.gain;
save_var.tx.gain_db = tx.gain_db;
save("results\rx_pattern.mat", "save_var");
