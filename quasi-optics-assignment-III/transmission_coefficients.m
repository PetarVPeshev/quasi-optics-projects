close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
medium.er = 11.9;
Ntheta = 1200;

%% DEPENDENT PARAMETERS
theta_crit = asin(1 / sqrt(medium.er));

%% SPHERICAL COORDINATE GRID
theta = linspace(eps, theta_crit, Ntheta);

%% TRANSMISSION ANGLE
theta_t = asin( sqrt(medium.er) .* sin(theta) );

%% CALCULATE TRANSMISSION COEFFICIENTS
[par_coeff, per_coeff] = transm_coeff(theta, theta_t, medium.er, 1);

%% CALCULATE TRANSMITTED POWER
[TE, TM] = surf_transm_power(par_coeff, per_coeff, theta, ...
    theta_t, medium.er, 1);

%% PLOT TRANSMISSION COEFFICIENTS
figure('Position', [250 250 700 400]);
plot(theta * 180 / pi, real(TE), 'LineWidth', 2.0, 'DisplayName', 'TE');
hold on;
plot(theta * 180 / pi, real(TM), 'LineWidth', 2.0, 'DisplayName', 'TM');
grid on;
xlim([0 20]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta_{i} / deg');
ylabel('P_{t} / P_{i}');
title('TE & TM Transmitted Power Ratio');
saveas(gcf, 'figures\transmitted_power.fig');
