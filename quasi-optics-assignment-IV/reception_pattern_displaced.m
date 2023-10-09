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

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;
feed.D = [4, 2, 1, 4] * wave.wavelength;
feed.dx = [4, 2, 1, 0] * wave.wavelength;
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

%% WAVE VECTOR
[k_comp, ~] = wave_vector(medium.er, wave.k, sph_grid);

%% RECEPTION GO FIELD
[rx.Ego, rx.Vgo, rx.sph_grid] = go_field(reflector.f, reflector.D, ...
    wave.E0, wave.k, sph_grid, sph_grid_i);

tx.Va = NaN( [size(sph_grid, 1, 2), 3, length(feed.D)] );
tx.rad_power = NaN(1, length(feed.D));
rx.Prx = NaN( [size(sph_grid_i, 1, 2), length(feed.D)] );
rx.Prx_norm = NaN( size(rx.Prx, 1, 2, 3) );
rx.eta = NaN( [size(sph_grid_i, 1, 2), length(feed.D)] );
rx.eta_max = NaN(1, length(feed.D));
for idx = 1 : 1 : length(feed.D)
    %% TRANSMISSION A FIELD
    [tx.Va(:, :, :, idx), tx.rad_power(idx)] = ...
        feed_va_field(feed.D(idx) / 2, wave.k, sph_grid, feed.dx(idx), ...
        k_comp(:, :, 1));
    
    %% RECEPTION POWER AND EFFICIENCY
    [rx.Prx(:, :, idx), rx.eta(:, :, idx)] = rx_power(rx.Vgo, ...
        tx.Va(:, :, :, idx), tx.rad_power(idx), wave.E0, ...
        reflector.D / 2, sph_grid, rx.sph_grid);
    rx.Prx_norm(:, :, idx) = norm_magnitude(rx.Prx(:, :, idx), 'dB');
    rx.eta_max(idx) = max(rx.eta(:, :, idx), [], 'all');
end

%% INCIDENT UV COORDINATES
uv_grid_i = uv_repr(sph_grid_i);

%% PLOT RECEPTION PATTERN IN UV REPRESENTATION
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(uv_grid_i(:, :, 1), uv_grid_i(:, :, 2), rx.Prx_norm(:, :, 1), ...
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
surface(uv_grid_i(:, :, 1), uv_grid_i(:, :, 2), rx.Prx_norm(:, :, 1), ...
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
    ' GHz, D = 100\lambda, D_{feed} = 4\lambda, f_{#} = 2,' ...
    ' and d_{x} = 4\lambda'], 'FontSize', 17, 'FontWeight', 'bold');
saveas(gcf, 'figures\rx_pattern_dx4.fig');

%% PLOT E PLANE RECEPTION AND TRANSMISSION PATTERNS
theta_plot = NaN(1, length(theta_i) * 2);
theta_plot(1 : length(theta_i)) = - fliplr(theta_i) * 180 / pi;
theta_plot(length(theta_i) + 1 : end) = theta_i * 180 / pi;
figure('Position', [250 250 725 400]);
plane_idx_1 = find(round(phi * 180 / pi, 0) == 0, 1);
plane_idx_2 = find(floor(phi * 180 / pi) == 180, 1);
plane_field = NaN(1, length(theta_i) * 2);
plane_field(1 : length(theta_i)) = fliplr(rx.Prx_norm(plane_idx_2, :, 1));
plane_field(length(theta_i) + 1 : end) = rx.Prx_norm(plane_idx_1, :, 1);
plot(theta_plot, plane_field, 'LineWidth', 2.0, ...
    'DisplayName', 'P_{rx}, D_{feed} = 4\lambda and d_{x} = 4\lambda');
hold on;
plane_field = NaN(1, length(theta_i) * 2);
plane_field(1 : length(theta_i)) = fliplr(rx.Prx_norm(plane_idx_2, :, 4));
plane_field(length(theta_i) + 1 : end) = rx.Prx_norm(plane_idx_1, :, 4);
plot(theta_plot, plane_field, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'P_{rx}, D_{feed} = 4\lambda and d_{x} = 0');
grid on;
ylim([-60 0]);
xlim([-3 3]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|P_{rx}| / dBs');
title(['RX Pattern @ E-plane, f = ' num2str(wave.f * 1e-9) ...
    ' GHz, D = 100\lambda, and f_{#} = 2']);
saveas(gcf, 'figures\rx_pattern_dx4_and_dx0.fig');

%% PLOT E PLANE RECEPTION AND TRANSMISSION PATTERNS FOR ALL DISPLACEMENTS
theta_plot = NaN(1, length(theta_i) * 2);
theta_plot(1 : length(theta_i)) = - fliplr(theta_i) * 180 / pi;
theta_plot(length(theta_i) + 1 : end) = theta_i * 180 / pi;
figure('Position', [250 250 725 400]);
plane_idx_1 = find(round(phi * 180 / pi, 0) == 0, 1);
plane_idx_2 = find(floor(phi * 180 / pi) == 180, 1);
for idx = 1 : 1 : length(feed.D)
    plane_field = NaN(1, length(theta_i) * 2);
    plane_field(1 : length(theta_i)) = ...
        fliplr(rx.Prx_norm(plane_idx_2, :, idx));
    plane_field(length(theta_i) + 1 : end) = ...
        rx.Prx_norm(plane_idx_1, :, idx);
    plot(theta_plot, plane_field, 'LineWidth', 2.0, ...
        'DisplayName', ['P_{rx}, D_{feed} = ' ...
        num2str(feed.D(idx) / wave.wavelength) '\lambda and d_{x} = ' ...
        num2str(feed.dx(idx) / wave.wavelength) '\lambda']);
    hold on;
end
hold off;
grid on;
ylim([-60 0]);
xlim([-3 3]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|P_{rx}| / dBs');
title(['RX Pattern @ E-plane, f = ' num2str(wave.f * 1e-9) ...
    ' GHz, D = 100\lambda, and f_{#} = 2']);
saveas(gcf, 'figures\rx_pattern_dx4_dx2_dx1_and_dx0.fig');

%% SAVE RX AND TX STRUCTS
save_var.eta = rx.eta;
save_var.eta_max = rx.eta_max;
save("results\rx_pattern_displaced.mat", "save_var");
