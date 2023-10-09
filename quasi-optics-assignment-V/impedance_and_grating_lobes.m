close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
wave.f = 10e9;
medium.er = 1;
dx = [15e-3 20e-3];
dy = [15e-3 20e-3];
dipole.w = 1e-3;
dipole.l = 14e-3;
Ntheta = 1200;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = [0 45 90 180 225 270] * pi / 180;
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, ~] = wave_vector(medium.er, wave.k0, sph_grid);

Zin = NaN( [size(sph_grid, 1, 2), 2] );
lobes = NaN( [9, 2001, 2, 2] );
for idx = 1 : 1 : length(dx)
    %% FLOQUET MODES
    [~, ~, modes] = floquet_modes(wave.k0, medium.er, ...
        dx(idx), dy(idx), sph_grid);

    %% ARRAY IMPEDANCE
    Zin(:, :, idx) = array_impedance(medium.er, dipole.l, dipole.w, ...
        dx(idx), dy(idx), wave.k0, k_comp);

    %% GRATING LOBES
    lobes(:, :, :, idx) = grating_lobes(wave.k0, dx(idx), dy(idx), modes);
end

%% PLOT E, D, AND H PLANE IMPEDANCE
theta_length = length(theta);
theta_plot = NaN(1, theta_length * 2);
theta_plot(1 : theta_length) = - fliplr(theta) * 180 / pi;
theta_plot(theta_length + 1 : end) = theta * 180 / pi;
color_styles = ["#0072BD", "#D95319"];
plane = ['E', 'D', 'H'];
y_limits = [-40 80; -50 500; -30 500];
for idx = 1 : 1 : 3
    figure('Position', [250 250 800 400]);
    plane_zin = NaN(1, theta_length * 2);
    plane_zin(1 : theta_length) = fliplr(Zin(idx + 3, :, 1));
    plane_zin(theta_length + 1 : end) = Zin(idx, :, 1);
    plot(theta_plot, real(plane_zin), 'LineWidth', 2.0, ...
        'Color', "#0072BD", 'DisplayName', '\Re(Z_{in})');
    hold on;
    plot(theta_plot, imag(plane_zin), 'LineWidth', 2.0, ...
        'Color', "#D95319", 'DisplayName', '\Im(Z_{in})');
    hold off;
    grid on;
    legend show;
    legend('location', 'bestoutside');
    xticks(-90 : 15 : 90);
    xlim([-90 90]);
    ylim(y_limits(idx, :));
    xlabel('\theta / deg');
    ylabel('Z_{in} / \Omega');
    title(['Array Impedance @ ' plane(idx) '-plane, f = ' ...
        num2str(wave.f * 1e-9) ' GHz, L = ' num2str(dipole.l * 1e3) ...
        ' mm, W = ' num2str(dipole.w * 1e3) ' mm, d_{x} = ' ...
        num2str(dx(1) * 1e3) ' mm, and d_{y} = ' num2str(dy(1) * 1e3) ...
        ' mm']);
    saveas(gcf, ['figures\Zin_d_' num2str(dx(1) * 1e3) 'mm_' plane(idx) ...
        '_plane.fig']);
end

y_limits = [-20 40; -50 500; -80 500];
for idx = 1 : 1 : 3
    figure('Position', [250 250 800 400]);
    plane_zin = NaN(1, theta_length * 2);
    plane_zin(1 : theta_length) = fliplr(Zin(idx + 3, :, 2));
    plane_zin(theta_length + 1 : end) = Zin(idx, :, 2);
    plot(theta_plot, real(plane_zin), 'LineWidth', 2.0, ...
        'Color', "#0072BD", 'DisplayName', '\Re(Z_{in})');
    hold on;
    plot(theta_plot, imag(plane_zin), 'LineWidth', 2.0, ...
        'Color', "#D95319", 'DisplayName', '\Im(Z_{in})');
    hold off;
    grid on;
    legend show;
    legend('location', 'bestoutside');
    xticks(-90 : 15 : 90);
    xlim([-90 90]);
    ylim(y_limits(idx, :));
    xlabel('\theta / deg');
    ylabel('Z_{in} / \Omega');
    title(['Array Impedance @ ' plane(idx) '-plane, f = ' ...
        num2str(wave.f * 1e-9) ' GHz, L = ' num2str(dipole.l * 1e3) ...
        ' mm, W = ' num2str(dipole.w * 1e3) ' mm, d_{x} = ' ...
        num2str(dx(2) * 1e3) ' mm, and d_{y} = ' num2str(dy(2) * 1e3) ...
        ' mm']);
    saveas(gcf, ['figures\Zin_d_' num2str(dx(2) * 1e3) 'mm_' plane(idx) ...
        '_plane.fig']);
end

%% PLOT GRATING LOBES
planes_theta = linspace(- pi / 2, pi / 2, 2001);
planes_k = NaN( [3, size(planes_theta)] );
planes_k(1, :) = wave.k0 * sin(planes_theta);
planes_k(2, :) = 0;
planes_k(3, 1) = wave.k0 * sin(45 * pi / 180);
planes_k(3, 2) = wave.k0 * sin(225 * pi / 180);
for idx = 1 : 1 : 2
    figure('Position', [250 250 600 400]);
    plot(lobes(:, :, 1, idx)', lobes(:, :, 2, idx)', ...
        'Color', [0 0 0], 'LineWidth', 2.0);
    hold on;
    patch(lobes(5, :, 1, idx), lobes(5, :, 2, idx), 'black', ...
        'FaceAlpha', 0.1);
    hold on;
    plot(planes_k(1, :), planes_k(2, :), ...
        'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2.0);
    hold on;
    plot(planes_k(2, :), planes_k(1, :), ...
        'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0);
    hold on;
    plot(planes_k(3, :), planes_k(3, :), ...
        'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0);
    grid on;
    legend show;
    legend('lobes', '', '', '', '', '', '', '', '', 'visible region', ...
        'E-plane', 'H-plane', 'D-plane', 'location', 'bestoutside');
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    xlabel('k_{mx}');
    ylabel('k_{my}');
    title(['Grating Lobes @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
        'd_{x} = ' num2str(dx(idx) * 1e3) ' mm, and d_{y} = ' ...
        num2str(dy(idx) * 1e3) ' mm']);
    saveas(gcf, ['figures\grating_lobes_' num2str(dx(idx) * 1e3) ...
        'mm.fig']);
end
