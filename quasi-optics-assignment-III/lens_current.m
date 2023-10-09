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
Nrho = 1200;
Nphi = 1200;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
lens.D = 6 * wave.wavelength;

%% CYLINDRICAL COORDINATE GRID
rho = linspace(eps, lens.D / 2, Nrho);
phi = linspace(eps, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);

%% APERTURE CURRENT
[Ja, lens.param, cart_grid] = lens_aperture(lens.D, wave.n, lens.er, ...
    medium.er, cyl_grid, lens.theta_max);

%% TRANSMITTED POWER
[par_coeff, per_coeff] = transm_coeff(lens.param.theta_i, ...
    lens.param.theta_t, medium.er, lens.er);
[TE, TM] = surf_transm_power(par_coeff, per_coeff, lens.param.theta_i, ...
    lens.param.theta_t, medium.er, lens.er);

%% NORMALIZED APERTURE CURRENT
Ja_max_magn = 20 * log10( max( max( max( abs(Ja) ) ) ) );
Ja_norm = NaN( [size(Ja, 1, 2, 3)] );
Ja_norm(:, :, 1) = 20 * log10( abs(Ja(:, :, 1)) ) - Ja_max_magn;
Ja_norm(:, :, 2) = 20 * log10( abs(Ja(:, :, 2)) ) - Ja_max_magn;
Ja_norm(:, :, 3) = 20 * log10( abs(Ja(:, :, 3)) );

%% PLOT TRANSMISSION COEFFICIENTS
figure('Position', [250 250 700 400]);
plot(lens.param.sph_grid(1, :, 2) * 180 / pi, real(TE(1, :)), ...
    'LineWidth', 2.0, 'DisplayName', 'TE');
hold on;
plot(lens.param.sph_grid(1, :, 2) * 180 / pi, real(TM(1, :)), ...
    'LineWidth', 2.0, 'DisplayName', 'TM');
grid on;
xticks(0 : 5 : 40);
yticks(0 : 0.2 : 1);
xlim([0 40]);
ylim([0 1]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta_{i} / deg');
ylabel('P_{t} / P_{i}');
title(['Lens Transmitted Power Ratio @ \phi = 0 deg, D = ' ...
    num2str(round(lens.D * 1e2, 2)) ' cm, \epsilon_{r} = ' ...
    num2str(lens.er) ', \theta_{0} = ' ...
    num2str(round(lens.theta_max * 180 / pi, 0)) ' deg']);
saveas(gcf, 'figures\lens_transmitted_power_aperture_current.fig');

%% PLOT ELECTRIC CURRENT DENSITY
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(cart_grid(:, :, 1) * 1e3, cart_grid(:, :, 2) * 1e3, ...
    Ja_norm(:, :, 1), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xticks(-5 : 2.5 : 5);
yticks(-5 : 2.5 : 5);
caxis([-10 0]);
view(0, 90);
xlabel('x / mm');
ylabel('y / mm');
zlabel('|J_{x}| / dB');
title('|J_{x}| / dB');
subplot(1, 2, 2);
surface(cart_grid(:, :, 1) * 1e3, cart_grid(:, :, 2) * 1e3, ...
    Ja_norm(:, :, 2), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
xticks(-5 : 2.5 : 5);
yticks(-5 : 2.5 : 5);
caxis([-45 -25]);
view(0, 90);
xlabel('x / mm');
ylabel('y / mm');
zlabel('|J_{y}| / dB');
title('|J_{y}| / dB');
sgtitle(['Aperture Electric Current Density @ f = ' ...
    num2str(wave.f * 1e-9) ' GHz, n = ' num2str(wave.n) ', D = ' ...
    num2str(round(lens.D * 1e2, 2)) ' cm, \epsilon_{r} = ' ...
    num2str(lens.er) ', \theta_{0} = ' ...
    num2str(round(lens.theta_max * 180 / pi, 0)) ' deg'], 'FontSize', ...
    17, 'FontWeight', 'bold');
saveas(gcf, 'figures\aperture_current.fig');
