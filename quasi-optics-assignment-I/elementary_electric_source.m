close all;
clear;
clc;

addpath('../../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
medium.er = 1;
wave.f = 28e9;
Nkx = 300;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;

%% WAVE VECTOR COMPONENTS
k_comp = zeros(1, Nkx, 3);
k_comp(:, :, 1) = wave.k0 * linspace(eps, 3, Nkx);
k_comp(:, :, 3) =  - 1j * sqrt( - wave.k0 .^ 2 + k_comp(:, :, 1) .^ 2 ...
    + k_comp(:, :, 2) .^ 2 );

%% CALCULATE EJ DYADIC SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(medium.er, wave.k0, k_comp, 'E', 'J');

%% PLOT X AND Y-COMPONENTS OF THE SPECTRAL GREEN'S FUNCTION
figure('Position', [250 250 700 400]);
plot(k_comp(:, :, 1) / wave.k0, real(SGFej(:, :, 1, 1)), ...
    'LineWidth', 3.0, 'DisplayName', '\Re\{G_{xx}^{EJ}\}');
hold on;
plot(k_comp(:, :, 1) / wave.k0, imag(SGFej(:, :, 1, 1)), ...
    'LineWidth', 3.0, 'DisplayName', '\Im\{G_{xx}^{EJ}\}');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('k_{x} / k_{0}');
ylabel('G_{xx}^{EJ}');
title("Spectral Green's Function XX-Component");

figure('Position', [250 250 700 400]);
plot(k_comp(:, :, 1) / wave.k0, real(SGFej(:, :, 2, 2)), ...
    'LineWidth', 3.0, 'DisplayName', '\Re\{G_{yy}^{EJ}\}');
hold on;
plot(k_comp(:, :, 1) / wave.k0, imag(SGFej(:, :, 2, 2)), ...
    'LineWidth', 3.0, 'DisplayName', '\Im\{G_{yy}^{EJ}\}');
grid on;
ylim([-1000 0]);
legend show;
legend('location', 'bestoutside');
xlabel('k_{x} / k_{0}');
ylabel('G_{yy}^{EJ}');
title("Spectral Green's Function YY-Component");
