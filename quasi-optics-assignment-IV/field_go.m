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
theta_i = 0;
phi_i = 0;
Nphi = 1200;
Ntheta = 300;
R = 1;
planes = [0, 90];

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;
feed.D = 4 * wave.wavelength;
reflector.D = 100 * wave.wavelength;
reflector.f = reflector.D * 2;

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - eps, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% INCIDENT ANGLE COORDINATE GRID
sph_grid_i(:, :, 1) = theta_i;
sph_grid_i(:, :, 2) = phi_i;

%% TRANSMITTED A FIELD
[tx.Va, tx.rad_power] = feed_va_field(feed.D / 2, wave.k, sph_grid);
tx.Va_total = total_field(tx.Va);

%% RECEIVED GO FIELD
[rx.Ego, rx.Vgo, rx.sph_grid] = go_field(reflector.f, reflector.D, ...
    wave.E0, wave.k, sph_grid, sph_grid_i);
rx.Vgo_total = total_field(rx.Vgo);

%% RECEIVED POWER AND EFFICIENCY
[rx.Prx, rx.eta] = rx_power(rx.Vgo, tx.Va, tx.rad_power, wave.E0, ...
    reflector.D / 2, sph_grid, rx.sph_grid);

%% PLOT GO AND A FIELD
theta_length = length( rx.sph_grid(1, :, 1) );
theta_plot = NaN(1, theta_length * 2);
theta_plot(1 : theta_length) = - fliplr(rx.sph_grid(1, :, 1)) * 180 / pi;
theta_plot(theta_length + 1 : end) = rx.sph_grid(1, :, 1) * 180 / pi;
line_styles = ["-", "--"];
figure('Position', [250 250 800 400]);
for idx = 1 : 1 : length(planes)
    plane_idx = find(round(phi * 180 / pi, 0) == planes(idx));
    plane_field = NaN(1, theta_length * 2);
    plane_field(1 : theta_length) = fliplr(rx.Vgo_total(plane_idx(1), :));
    plane_field(theta_length + 1 : end) = rx.Vgo_total(plane_idx(1), :);
    plot(theta_plot, norm_magnitude(plane_field, 'dB'), ...
        [line_styles(idx)], 'LineWidth', 3.0, 'DisplayName', ...
        ['V^{GO}_{rx}, \phi = ' num2str( round(planes(idx), 0) ) ...
        ' and ' num2str( round(planes(idx), 0) + 180 ) ' deg']);
    hold on;
end
theta_plot = NaN(1, length(theta) * 2);
theta_plot(1 : length(theta)) = - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
for idx = 1 : 1 : length(planes)
    plane_idx = find(round(phi * 180 / pi, 0) == planes(idx));
    plane_field = NaN(1, length(theta) * 2);
    plane_field(1 : length(theta)) = fliplr(tx.Va_total(plane_idx(1),:));
    plane_field(length(theta) + 1 : end) = tx.Va_total(plane_idx(1), :);
    plot(theta_plot, norm_magnitude(plane_field, 'dB'), ...
        [line_styles(idx)], 'LineWidth', 3.0, 'DisplayName', ...
        ['V^{A}_{tx}, \phi = ' num2str( round(planes(idx), 0) ) ' and ' ...
        num2str( round(planes(idx), 0) + 180 ) ' deg']);
end
hold off;
grid on;
legend show;
legend('location', 'bestoutside');
xticks(-15 : 5 : 15);
yticks(-5 : 1 : 0);
ylim([-5 0]);
xlim([- max(rx.sph_grid(1, :, 1)) max(rx.sph_grid(1, :, 1))] * 180 / pi);
xlabel('\theta / deg');
ylabel('|E| / dB');
title(['TX and RX Fields @ E^{PW}_{i} = ' num2str(wave.E0) ' V /m, ' ...
    '\theta_{i} = ' num2str(theta_i * 180 / pi) ' deg, and ' ...
    '\phi_{i} = ' num2str(phi_i * 180 / pi) ' deg']);
saveas(gcf, 'figures\rx_tx_fields.fig');

%% SAVE RECEIVED POWER AND APERTURE EFFICIENCY
save_var.Prx = rx.Prx;
save_var.eta = rx.eta;
save("results\rx_power_and_efficiency.mat", "save_var");
