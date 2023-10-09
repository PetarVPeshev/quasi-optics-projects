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
lens.er = 11.9;
Nrho = 1200;
Nphi = 1200;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
lens.D = 6 * wave.wavelength;

%% CYLINDRICAL COORDINATE GRID
rho = linspace(eps, lens.D / 2, Nrho);
phi = linspace(eps, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);

%% LENS APERTURE PARAMETERS
[lens.param] = lens_parameters(lens.D, lens.er, cyl_grid);

%% CALCULATE TRANSMISSION COEFFICIENTS
[par_coeff, per_coeff] = transm_coeff(lens.param.theta_i, ...
    lens.param.theta_t, lens.er, 1);

%% CALCULATE TRANSMITTED POWER
[TE_power, TM_power] = surf_transm_power(par_coeff, per_coeff, ...
    lens.param.theta_i, lens.param.theta_t, lens.er, 1);

%% PLOT TRANSMISSION COEFFICIENTS
figure('Position', [250 250 700 400]);
plot(lens.param.sph_grid(1, :, 2) * 180 / pi, real(TE_power(1, :)), ...
    'LineWidth', 2.0, 'DisplayName', 'TE');
hold on;
plot(lens.param.sph_grid(1, :, 2) * 180 / pi, real(TM_power(1, :)), ...
    'LineWidth', 2.0, 'DisplayName', 'TM');
grid on;
xticks(0 : 5 : 80);
xlim([0 80]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('P_{t} / P_{i}');
title(['Lens TE & TM Transmitted Power @ \phi = 0 deg, D = ' ...
    num2str(round(lens.D * 1e2, 2)) ' cm, \epsilon_{r} = ' ...
    num2str(lens.er)]);
saveas(gcf, 'figures\lens_transmitted_power.fig');
