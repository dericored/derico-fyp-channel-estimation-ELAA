clc; 
clear all
close all 

sparsity=20;
%%% system parameters
N = 512; % number of beams (transmit antennas)
L = 12; % number of all paths
gamma=0.5; 
Lf = L*gamma; % number of paths for far-field 
Ln = L*(1-gamma); % number of paths for near-field
M = 256; % number of pilot overhead

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space

% the far-field angle-domain DFT matrix
Uf = (1/sqrt(N))*exp(-1i*pi*[0:N-1]'*[-(N-1)/2:1:(N/2)]*(2/N));

% the near-field polar-domain transform matrix [4]
Rmin=10;
Rmax=80;
eta = 2.5; 
[Un, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);

% generate the far-field path components
hf = zeros(N,1);
alpha_f = (randn(Lf,1) + 1j*randn(Lf,1))/sqrt(2); % the gain 
sector = 3*pi/3; 
theta_f = rand(Lf,1) * sector - sector/2; % the angle [-pi/2, pi/2]
for l = 1:Lf
    af = far_field_manifold(N,theta_f(l));
    hf = hf + alpha_f(l)*af;
end

% generate the near-field path components
hn = zeros(N,1);
r_n = rand(Ln,1) * (Rmax - Rmin) + Rmin; % the distance 
sector = 3*pi/3; 
theta_n = rand(Ln,1) * sector - sector/2; % the angle [-pi/2, pi/2]
alpha_n = (randn(Ln,1) + 1j*randn(Ln,1))/sqrt(2); % the gain 
for l = 1:Ln
   an = near_field_manifold( N, d, fc, r_n(l), theta_n(l));
   hn = hn + alpha_n(l) * an;
end 

% hn in polar and angle domain
[hn_polar, ~] = OMP(hn, Un, sparsity*Ln);
hn_polar_power = abs(hn_polar);

hn_angle = Uf\hn;
hn_angle_power = abs(hn_angle);

% hf in polar and angle domain
[hf_polar, ~] = OMP(hf, Un, sparsity*Lf);
hf_polar_power = abs(hf_polar);

hf_angle = pinv(Uf)*hf;
hf_angle_power = abs(hf_angle);

% hf antenna vs. angle vs. polar domain
figure(1);
subplot(3,1,1);
stem(abs(hf), 'filled', 'MarkerSize', 1.5, 'color', '#27346F');
title("Non-sparse hf");
xlim([0 512]);
grid on;
subplot(3,1,2);
stem(abs(hf_angle), 'filled', 'color', [0.6350 0.0780 0.1840], 'MarkerSize', 1.5);
title("Angle-domain hf");
xlim([0 512]);
grid on;
subplot(3,1,3);
stem(abs(flip(hf_polar)), 'filled', 'color', [0.9290 0.6940 0.1250], 'MarkerSize', 1.5);
title("Polar-domain hf");
xlim([0 length(hf_polar)]);
xlabel("Vector element");
grid on;

% hn antenna vs. angle vs. polar domain
figure(2);
subplot(3,1,1);
stem(abs(hn), 'filled', 'MarkerSize', 1.5, 'color', '#27346F');
title("Non-sparse hn");
xlim([0 512]);
grid on;
subplot(3,1,2);
stem(abs(hn_angle), 'filled', 'color', [0.6350 0.0780 0.1840], 'MarkerSize', 1.5);
title("Angle-domain hn");
xlim([0 length(hn_angle)]);
grid on;
subplot(3,1,3);
stem(abs(flip(hn_polar)), 'filled', 'color', [0.9290 0.6940 0.1250], 'MarkerSize', 1.5);
title("Polar-domain hn");
xlim([0 length(hn_polar)]);
xlabel("Vector element");
grid on;