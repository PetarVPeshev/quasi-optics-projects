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
Nphi = 2000;
Ntheta = 2000;
R = 1;
planes = [0, 90];

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;
feed.D = 4 * wave.wavelength;

%% SPHIRICAL COORDINATES GRID
theta = linspace(eps, pi / 2 - eps, Ntheta);
phi = linspace(eps, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, ~] = wave_vector(medium.er, wave.k, sph_grid);

%% DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k, k_comp, 'E', 'J');

%% CIRCULAR FEED CURRENT DENSITY FOURIER TRANSFORM
Jft = ft_current(wave.k, feed.D / 2, sph_grid(:, :, 1), 'circular', 'y');
Jft_total = total_field(Jft);

%% ELECTRIC FAR-FIELD OF CIRCULAR FEED
Efar = farfield(wave.k, R, sph_grid, k_comp(:, :, 3), SGFej, Jft);
Efar_total = total_field(Efar);

%% UV COORDINATES
uv_grid = uv_repr(sph_grid);

%% PLOT FEED ELECTRIC CURRENT DENSITY
theta_plot = NaN(1, length(theta) * 2);
theta_plot(1 : length(theta)) = - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
line_styles = ["-"; "--"];
figure('Position', [250 250 725 400]);
for idx = 1 : 1 : length(planes)
    plane_idx = find(round(phi * 180 / pi, 1) == planes(idx), 1);
    plane_jft = NaN(1, length(theta) * 2);
    plane_jft(1 : length(theta)) = fliplr(Jft_total(plane_idx, :));
    plane_jft(length(theta) + 1 : end) = Jft_total(plane_idx, :);
    plot(theta_plot, norm_magnitude(plane_jft, 'dB'), ...
        [line_styles(idx)], 'LineWidth', 3.0, ...
        'DisplayName', ['|J_{FT}|, \phi = ' num2str(planes(idx)) ...
        ' deg and ' num2str(planes(idx) + 180) ' deg']);
    hold on;
end
hold off;
grid on;
ylim([-40 0]);
xlim([-90 90]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|J_{FT}| / dB');
title(['Feed J_{FT} @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz and D_{feed} = 4\lambda']);
saveas(gcf, 'figures\feed_J_FT.fig');

%% PLOT FEED ELECTRIC FAR-FIELD IN UV REPRESENTATION
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
caxis([-40 0]);
view(0, 90);
xlabel('U');
ylabel('V');
zlabel('|E_{feed}| / dB');
subplot(1, 2, 2);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
zlim([-40 0]);
caxis([-40 0]);
view(-37.5, 30);
xlabel('U');
ylabel('V');
zlabel('|E_{feed}| / dB');
sgtitle(['Feed E Far-Field @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz and D_{feed} = 4\lambda'], 'FontSize', 17, 'FontWeight', 'bold');
saveas(gcf, 'figures\feed_pattern.fig');

%% PLOT FEED ELECTRIC FAR-FIELD AT PHI 0, 90 DEGREES
theta_plot = NaN(1, length(theta) * 2);
theta_plot(1 : length(theta)) = - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
line_styles = ["-", "--"];
figure('Position', [250 250 725 400]);
for idx = 1 : 1 : length(planes)
    plane_idx = find(round(phi * 180 / pi, 1) == planes(idx), 1);
    plane_field = NaN(1, length(theta) * 2);
    plane_field(1 : length(theta)) = fliplr(Efar_total(plane_idx, :));
    plane_field(length(theta) + 1 : end) = Efar_total(plane_idx, :);
    plot(theta_plot, norm_magnitude(plane_field, 'dB'), ...
        [line_styles(idx)], 'LineWidth', 3.0, ...
        'DisplayName', ['E, \phi = ' num2str(planes(idx)) ' deg and ' ...
        num2str(planes(idx) + 180) ' deg']);
    hold on;
end
hold off;
grid on;
ylim([-40 0]);
xlim([-90 90]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E_{feed}| / dB');
title(['Feed E Far-Field @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz and D_{feed} = 4\lambda']);
saveas(gcf, 'figures\feed_E_planes.fig');
