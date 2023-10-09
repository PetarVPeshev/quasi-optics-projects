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
medium.er = 1;
lens.er = 11.9;
lens.theta_max = 40 * pi / 180;
Nrho = 200;
Nphi = 200;
Ntheta = 200;
R = 1;
planes = [0, 90];

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k = 2 * pi / wave.wavelength;
lens.D = 6 * wave.wavelength;

%% CYLINDRICAL COORDINATE GRID
rho = linspace(eps, lens.D / 2, Nrho);
phi = linspace(eps, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);

%% SPHERICAL COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.001, Ntheta);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, ~] = wave_vector(medium.er, wave.k, sph_grid);

%% DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k, k_comp, 'E', 'J');

%% APERTURE CURRENT
[~, lens.param, cart_grid, Jft] = lens_aperture(lens.D, wave.n, ...
    lens.er, medium.er, cyl_grid, lens.theta_max, 'FT', k_comp);

%% APERTURE ELECTRIC FAR-FIELD
lens.Efar = farfield(wave.k, R, sph_grid, k_comp(:, :, 3), SGFej, Jft);
lens.Efar_total = total_field(lens.Efar);

%% APERTURE RADIATION INTENSITY, RADIATED POWER, AND DIRECTIVITY
[lens.dir, lens.rad_intensity, lens.rad_power] = directivity(medium.er, ...
    lens.Efar, sph_grid, R);

%% FEED RADIATED POWER
[~, feed.rad_power] = lens_feed(wave.n, R, lens.er, sph_grid, wave.f);

%% EFFICIENCY AND GAIN
lens.figures.eff = lens.rad_power / feed.rad_power;
lens.figures.G = lens.dir * lens.figures.eff;
lens.figures.G_db = 10 * log10(lens.figures.G);
lens.figures.Gmax = max( max(lens.figures.G) );
lens.figures.Gmax_db = 10 * log10(lens.figures.Gmax);

%% UNIFORM APERTURE CURRENT DENSITY FOURIER TRANSFORM
uniform.Jft = ft_current(wave.k, lens.D / 2, sph_grid(:, :, 1), ...
    'circular', 'y');

%% UNIFORM APERTURE ELECTRIC FAR-FIELD
uniform.Efar = farfield(wave.k, R, sph_grid, k_comp(:, :, 3), ...
    SGFej, uniform.Jft);
uniform.Efar_total = total_field(uniform.Efar);

%% UV COORDINATES
uv_grid = uv_repr(sph_grid);

%% PLOT APERTURE ELECTRIC FAR-FIELD
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(lens.Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([-1 1]);
ylim([-1 1]);
caxis([-40 0]);
view(0, 90);
xticks(-1 : 1 : 1);
yticks(-1 : 0.5 : 1);
xlabel('U');
ylabel('V');
zlabel('|E| / dB');
subplot(1, 2, 2);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(lens.Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([-1 1]);
ylim([-1 1]);
caxis([-40 0]);
view(-37.5, 30);
xticks(-1 : 1 : 1);
yticks(-1 : 0.5 : 1);
xlabel('U');
ylabel('V');
zlabel('|E| / dB');
sgtitle(['E Far-Field @ f = ' num2str(wave.f * 1e-9) ' GHz, n = ' ...
    num2str(wave.n) ', D = ' num2str(round(lens.D * 1e2, 2)) ...
    ' cm, \epsilon_{r} = ' num2str(lens.er) ', \theta_{0} = ' ...
    num2str(round(lens.theta_max * 180 / pi, 0)) ' deg'], 'FontSize', ...
    17, 'FontWeight', 'bold');
saveas(gcf, 'figures\pattern.fig');

%% PLOT UNIFORM APERTURE ELECTRIC FAR-FIELD
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(uniform.Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([-1 1]);
ylim([-1 1]);
caxis([-40 0]);
view(0, 90);
xticks(-1 : 1 : 1);
yticks(-1 : 0.5 : 1);
xlabel('U');
ylabel('V');
zlabel('|E| / dB');
subplot(1, 2, 2);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(uniform.Efar_total, 'dB'), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xlim([-1 1]);
ylim([-1 1]);
caxis([-40 0]);
view(-37.5, 30);
xticks(-1 : 1 : 1);
yticks(-1 : 0.5 : 1);
xlabel('U');
ylabel('V');
zlabel('|E| / dB');
sgtitle(['Uniform Aperture E Far-Field @ f = ' num2str(wave.f * 1e-9) ...
    ' GHz, D = ' num2str(round(lens.D * 1e2, 2)) ' cm'], ...
    'FontSize', 17, 'FontWeight', 'bold');
saveas(gcf, 'figures\uniform_pattern.fig');

%% PLOT APERTURE AND UNIFORM APERTURE ELECTRIC FAR-FIELD AT PHI 0, 90 DEG
theta_length = length(theta);
theta_plot = NaN(1, theta_length * 2);
theta_plot(1 : theta_length) = - fliplr(theta) * 180 / pi;
theta_plot(theta_length + 1 : end) = theta * 180 / pi;
line_styles = ["-", "--"];
figure('Position', [250 250 725 400]);
for idx = 1 : 1 : length(planes)
    plane_idx = find(round(phi * 180 / pi, 0) == planes(idx), 1);
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
    plane_idx = find(round(phi * 180 / pi, 0) == planes(idx), 1);
    plane_field = NaN(1, theta_length * 2);
    plane_field(1 : theta_length) = fliplr(lens.Efar_total(plane_idx, :));
    plane_field(theta_length + 1 : end) = lens.Efar_total(plane_idx, :);
    plot(theta_plot, norm_magnitude(plane_field, 'dB'), ...
        [line_styles(idx)], 'LineWidth', 3.0, ...
        'DisplayName', ['\phi = ' num2str(planes(idx)) ' deg and ' ...
        num2str(planes(idx) + 180) ' deg']);
    hold on;
end
hold off;
grid on;
ylim([-40 0]);
xlim([-50 50]);
xticks(-50 : 10 : 50);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / dB');
title(['E Far-Field @ f = ' num2str(wave.f * 1e-9) ' GHz, n = ' ...
    num2str(wave.n) ', D = ' num2str(round(lens.D * 1e2, 2)) ...
    ' cm, \epsilon_{r} = ' num2str(lens.er) ', \theta_{0} = ' ...
    num2str(round(lens.theta_max * 180 / pi, 0)) ' deg']);
saveas(gcf, 'figures\E_planes.fig');
