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
feed.Nphi = 1000;
feed.Ntheta = 1000;
reflector.Nphi = 400;
reflector.Ntheta = 400;
planes = [0, 90];
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;
feed.D = 4 * wave.wavelength;
reflector.D = reflector.f / 2;

%% FEED COORDINATES GRID
feed.theta = linspace(eps, pi / 2 - 0.1 * pi / 180, feed.Ntheta);
feed.phi = linspace(eps, 2 * pi, feed.Nphi);
feed.sph_grid = meshgrid_comb(feed.theta, feed.phi);

%% REFLECTOR COORDINATES GRID
reflector.theta = linspace(eps, 4 * wave.wavelength / reflector.D, ...
    reflector.Ntheta);
reflector.phi = linspace(eps, 2 * pi, reflector.Nphi);
reflector.sph_grid = meshgrid_comb(reflector.theta, reflector.phi);

%% WAVE VECTOR COMPONENTS
[feed.k_comp, ~] = wave_vector(medium.er, wave.k, feed.sph_grid);
[reflector.k_comp, ~] = wave_vector(medium.er, wave.k, reflector.sph_grid);

%% DYADIC SPECTRAL GREEN'S FUNCTION
feed.SGFej = dyadic_sgf(medium.er, wave.k, feed.k_comp, 'E', 'J');
reflector.SGFej = dyadic_sgf(medium.er, wave.k, reflector.k_comp, ...
    'E', 'J');

%% FEED CURRENT DENSITY FOURIER TRANSFORM
feed.Jft = ft_current(wave.k, feed.D / 2, feed.sph_grid(:, :, 1), ...
    'circular', 'y');

%% FEED ELECTRIC FAR-FIELD
feed.Ecomp = farfield(wave.k, R, feed.sph_grid, feed.k_comp(:, :, 3), ...
    feed.SGFej, feed.Jft) * 2 * pi * R / exp(-1j * wave.k * R);

%% APERTURE CURRENT FOURIER TRANSFORM
[~, ~, reflector.Jft] = aperture_current(feed.Ecomp, wave.k, medium.er, ...
    feed.sph_grid, 'reflector', reflector.f, reflector.D, ...
    'Schelkunoff', 'FT', reflector.k_comp);

%% APERTURE ELECTRIC FAR-FIELD
reflector.Efar = farfield(wave.k, R, reflector.sph_grid, ...
    reflector.k_comp(:, :, 3), reflector.SGFej, reflector.Jft);
reflector.Efar_total = total_field(reflector.Efar);

%% UNIFORM APERTURE CURRENT DENSITY FOURIER TRANSFORM
uniform.Jft = ft_current(wave.k, reflector.D / 2, ...
    reflector.sph_grid(:, :, 1), 'circular', 'y');

%% UNIFORM APERTURE ELECTRIC FAR-FIELD
uniform.Efar = farfield(wave.k, R, reflector.sph_grid, ...
    reflector.k_comp(:, :, 3), reflector.SGFej, uniform.Jft);
uniform.Efar_total = total_field(uniform.Efar);

%% UV COORDINATES
reflector.uv_grid = uv_repr(reflector.sph_grid);

%% PLOT APERTURE ELECTRIC FAR-FIELD
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(reflector.uv_grid(:, :, 1), reflector.uv_grid(:, :, 2), ...
    norm_magnitude(reflector.Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([- max(reflector.theta) max(reflector.theta)]);
ylim([- max(reflector.theta) max(reflector.theta)]);
caxis([-40 0]);
view(0, 90);
xlabel('U');
ylabel('V');
zlabel('|E| / dB');
subplot(1, 2, 2);
surface(reflector.uv_grid(:, :, 1), reflector.uv_grid(:, :, 2), ...
    norm_magnitude(reflector.Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([- max(reflector.theta) max(reflector.theta)]);
ylim([- max(reflector.theta) max(reflector.theta)]);
caxis([-40 0]);
view(-37.5, 30);
xlabel('U');
ylabel('V');
zlabel('|E| / dB');
sgtitle(['E Far-Field @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    'D_{feed} = 4\lambda, F = ' num2str(reflector.f) ' m, and ' ...
    'D = ' num2str(reflector.D) ' m'], 'FontSize', 17, ...
    'FontWeight', 'bold');
saveas(gcf, 'figures\pattern.fig');

%% PLOT UNIFORM APERTURE ELECTRIC FAR-FIELD
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(reflector.uv_grid(:, :, 1), reflector.uv_grid(:, :, 2), ...
    norm_magnitude(uniform.Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([- max(reflector.theta) max(reflector.theta)]);
ylim([- max(reflector.theta) max(reflector.theta)]);
caxis([-40 0]);
view(0, 90);
xlabel('U');
ylabel('V');
zlabel('|E| / dB');
subplot(1, 2, 2);
surface(reflector.uv_grid(:, :, 1), reflector.uv_grid(:, :, 2), ...
    norm_magnitude(uniform.Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([- max(reflector.theta) max(reflector.theta)]);
ylim([- max(reflector.theta) max(reflector.theta)]);
caxis([-40 0]);
view(-37.5, 30);
xlabel('U');
ylabel('V');
zlabel('|E| / dB');
sgtitle(['Uniform Aperture E Far-Field @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz, F = ' num2str(reflector.f) ' m, and D = ' ...
    num2str(reflector.D) ' m'], 'FontSize', 17, 'FontWeight', 'bold');
saveas(gcf, 'figures\uniform_pattern.fig');

%% PLOT APERTURE AND UNIFORM APERTURE ELECTRIC FAR-FIELD AT PHI 0, 90 DEG
theta_length = length(reflector.theta);
theta_plot = NaN(1, theta_length * 2);
theta_plot(1 : theta_length) = - fliplr(reflector.theta) * 180 / pi;
theta_plot(theta_length + 1 : end) = reflector.theta * 180 / pi;
line_styles = ["-", "--"];
figure('Position', [250 250 725 400]);
for idx = 1 : 1 : length(planes)
    plane_idx = find(round(reflector.phi * 180 / pi, 0) == planes(idx), 1);
    plane_field = NaN(1, theta_length * 2);
    plane_field(1 : theta_length) = ...
        fliplr(uniform.Efar_total(plane_idx, :));
    plane_field(theta_length + 1 : end) = uniform.Efar_total(plane_idx, :);
    plot(theta_plot, norm_magnitude(plane_field, 'dB'), ...
        [line_styles(idx)], 'LineWidth', 3.0, ...
        'DisplayName', ['uniform, \phi = ' num2str(planes(idx)) ...
        ' deg and ' num2str(planes(idx) + 180) ' deg']);
    hold on;
end
for idx  = 1 : 1 : length(planes)
    plane_idx = find(round(reflector.phi * 180 / pi, 0) == planes(idx), 1);
    plane_field = NaN(1, theta_length * 2);
    plane_field(1 : theta_length) = ...
        fliplr(reflector.Efar_total(plane_idx, :));
    plane_field(theta_length + 1 : end) = ...
        reflector.Efar_total(plane_idx, :);
    plot(theta_plot, norm_magnitude(plane_field, 'dB'), ...
        [line_styles(idx)], 'LineWidth', 3.0, ...
        'DisplayName', ['\eta_{tap} = 0.83, \phi = ' ...
        num2str(planes(idx)) ' deg and ' num2str(planes(idx) + 180) ...
        ' deg']);
    hold on;
end
hold off;
grid on;
ylim([-40 0]);
xlim([min(theta_plot) max(theta_plot)]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / dB');
title(['E Far-Field @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    'D_{feed} = 4\lambda, F = ' num2str(reflector.f) ' m, and D = ' ...
    num2str(reflector.D) ' m']);
saveas(gcf, 'figures\E_planes.fig');
