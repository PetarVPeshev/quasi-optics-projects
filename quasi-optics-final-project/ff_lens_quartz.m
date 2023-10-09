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
zeta = 376.730313668;

%% PARAMETERS
% Wave parameters
wave.f = 70 * 1e9;
% Waveguide parameters
waveguide.E10 = 1;
waveguide.er = 1;
waveguide.a = 2.50 * 1e-3;
waveguide.b = 1.88 * 1e-3;
% Lens parameters
lens.er = 4;
lens.theta_max = 35 * pi / 180;
% Grids parameters
Nrho = 100;
Nphi = 100;
Ntheta = 100;
% Far-field parameters
R = 1;

%% DEPENDENT PARAMETERS
% Wave parameters
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
% Waveguide parameters
waveguide.kx1 = pi / waveguide.a;
waveguide.kz = sqrt(wave.k0 ^ 2 - waveguide.kx1 .^ 2);
waveguide.ZTE = zeta * wave.k0 ./ waveguide.kz;
waveguide.n = sqrt(waveguide.er);
% Lens parameters
lens.D = 10 * wave.wavelength;
lens.n = sqrt(lens.er);
lens.Z = zeta / lens.n;

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
waveguide_lens.TE_coef = 2 * lens.Z ./ (lens.Z + waveguide.ZTE);

%% CYLINDRICAL COORDINATE GRID INSIDE LENS
rho = linspace(eps, lens(1).D / 2, Nrho);
phi = linspace(0, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);

%% SPHERICAL COORDINATE GRID OUTSIDE LENS
theta = linspace(eps,  pi / 2 - 0.1 * pi / 180, Ntheta);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k0_comp, ~] = wave_vector(1, wave.k0, sph_grid);

%% SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(1, wave.k0, k0_comp, 'E', 'J');

%% LENS PARAMETERS
[lens, ~, lens_sph_grid, theta_inc, theta_tr] ...
    = lens_calculation(lens, cyl_grid, lens.theta_max);
    
%% RADIAL DISTANCE TO LENS INTERFACE AND SPHERICAL GRID INSIDE LENS
r = lens_sph_grid(:, :, 1);
lens_sph_grid = lens_sph_grid(:, :, 2 : 3);

%% WAVE VECTOR COMPONENTS
[k_comp, k] = wave_vector(lens.er, wave.k0, lens_sph_grid);

%% WAVEGUIDE FEED
[Ewg, Mwg] = waveguide_feed(waveguide, waveguide_lens.TE_coef, ...
    k, k_comp, r, lens_sph_grid, 'NeglectPhase');

%% APERTURE CURRENT
[J, lens, Jft] = lens_current(lens, Ewg, theta_inc, theta_tr, 1, ...
    cyl_grid, lens_sph_grid, 'FT', k0_comp);

%% RADIATED FAR-FIELD
E = farfield(wave.k0, R, sph_grid, k0_comp(:, :, 3), SGFej, Jft);
E_total = total_field(E);

%% APERTURE DIRECTIVITY, AND RADIATED POWER
[dir, ~, rad_power] = directivity(1, E, sph_grid, R);
dir_db = 10 * log10(dir);
dir_broadside_db = dir_db(1, 1);

%% EFFICIENCY
[k_comp_feed, k_feed] = wave_vector(lens.er, wave.k0, sph_grid);
[Ewg_rad, ~] = waveguide_feed(waveguide, waveguide_lens.TE_coef, ...
    k_feed, k_comp_feed, R, sph_grid);
[~, ~, rad_power_feed] = directivity(lens.er, Ewg_rad, sph_grid, R);
eta = rad_power / rad_power_feed;

%% CARTESIAN COORDINATES
cyl_coord = zeros( [size(cyl_grid, 1, 2), 3] );
cyl_coord(:, :, 1 : 2) = cyl_grid;
cart_coord = cyl2cart_cord(cyl_coord);

%% PLOT EQUIVALENT APERTURE ELECTRIC CURRENT DENSITY
J_max_magn = 20 * log10( max(abs(J), [], 'all') );
J_norm = NaN( [size(cart_coord, 1, 2), 3] );
J_norm(:, :, 1) = 20 * log10(abs(J(:, :, 1))) - J_max_magn;
J_norm(:, :, 2) = 20 * log10(abs(J(:, :, 2))) - J_max_magn;
J_norm(:, :, 3) = 20 * log10(abs(J(:, :, 3))) - J_max_magn;

figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(cart_coord(:, :, 1) * 1e3, cart_coord(:, :, 2) * 1e3, ...
    J_norm(:, :, 1), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
caxis([-50 0]);
view(0, 90);
xlabel('x / mm');
ylabel('y / mm');
zlabel('|J_{x}| / dB');
title('|J_{x}| / dB');
subplot(1, 2, 2);
surface(cart_coord(:, :, 1) * 1e3, cart_coord(:, :, 2) * 1e3, ...
    J_norm(:, :, 2), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
caxis([-50 0]);
view(0, 90);
xlabel('x / mm');
ylabel('y / mm');
zlabel('|J_{y}| / dB');
title('|J_{y}| / dB');
sgtitle(['Aperture J_{eq} @ quartz, f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    'D = 10\lambda, \theta_{max} = 35 deg, a = ' ...
    num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
    num2str(round(waveguide.b * 1e3, 2)) ' mm'], ...
    'FontSize', 17, 'FontWeight', 'bold');
saveas(gcf, 'figures\jeq_lens_quartz.fig');

%% UV COORDINATES
uv_grid = uv_repr(sph_grid);

%% PLOT FAR-FIELD
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(E_total, 'dB'), 'LineStyle', 'none');
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
    norm_magnitude(E_total, 'dB'), 'LineStyle', 'none');
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
sgtitle(['E^{FF} @ quartz, f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    'D = 10\lambda, \theta_{max} = 35 deg, a = ' ...
    num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
    num2str(round(waveguide.b * 1e3, 2)) ' mm'], ...
    'FontSize', 17, 'FontWeight', 'bold');
saveas(gcf, 'figures\ff_lens_quartz.fig');

%% SAVE WORKSPACE
lens_quartz.sph_grid = sph_grid;
lens_quartz.theta = theta;
lens_quartz.phi = phi;
lens_quartz.cart_coord = cart_coord;
lens_quartz.E = E;
lens_quartz.J = J;
lens_quartz.dir_broadside_db = dir_broadside_db;
lens_quartz.eta = eta;
save('results\lens_quartz.mat', 'lens_quartz');
