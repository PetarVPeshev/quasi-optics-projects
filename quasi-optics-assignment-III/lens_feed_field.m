close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
wave.f = 180e9;
wave.n = 4;
medium.er = 11.9;
Nphi = 1200;
Ntheta = 300;
R = 1;

%% SPHERICAL COORDINATE GRID
theta = linspace(eps, pi / 2 - eps, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% CALCULATE FEED ELECTRIC FAR-FIELD
[Efeed, Prad_feed] = lens_feed(wave.n, R, medium.er, sph_grid, wave.f);
Efeed_total = total_field(Efeed);

%% UV COORDINATES
uv_mesh = uv_repr(sph_grid);

%% PLOT ELECTRIC FAR-FIELD USING UV REPRESENTATION
Efeed_norm = norm_magnitude(Efeed_total);
figure();
surface(uv_mesh(:, :, 1), uv_mesh(:, :, 2), Efeed_norm, ...
    'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
caxis([-40 0]);
xlabel('U');
ylabel('V');
zlabel('|E_{feed}| / dB');
title(['Feed E Far-Field @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    '\epsilon_{r} = ' num2str(medium.er) ', n = ' num2str(wave.n) ]);
saveas(gcf, 'figures\lens_feed.fig');

plot_idx = find(round(phi * 180 / pi, 0) == 180, 1);
theta_plot = NaN(1, length(theta) * 2);
theta_plot(1 : length(theta)) = - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
plane_plot = NaN(1, length(theta) * 2);
plane_plot(1 : length(theta)) = fliplr(Efeed_norm(plot_idx, :));
plane_plot(length(theta) + 1 : end) = Efeed_norm(1, :);
figure('Position', [250 250 725 400]);
plot(theta_plot, plane_plot, 'LineWidth', 2.0, ...
    'DisplayName', 'E_{feed}, \phi = 0 deg and 180 deg');
grid on;
xticks(-90 : 15 : 90);
xlim([-90 90]);
ylim([-40 0]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E_{feed}| / dB');
title(['Feed E Far-Field @ E-plane, f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    '\epsilon_{r} = ' num2str(medium.er) ', n = ' num2str(wave.n) ]);
saveas(gcf, 'figures\lens_feed_E_plane.fig');
