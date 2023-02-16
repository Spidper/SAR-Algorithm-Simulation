clear;
c = 3e8;            % Speed of light
Rc = 20000;         % Range of image center
V = 150;            % Speed of radar
Tr = 2.5e-6;        % Pulse width
Kr = 20e12;         % Frequency modulation rate of LFM signal
f0 = 5.3e9;         % Center frequency
Fr = 60e6;          % Range sampling frequency
Fa = 100;           % Azimuth sampling frequency
Nr = 256;           % Range sampled points
Na = 256;           % Azimuth sampled points
Theta_sq = 10 * pi/180;  % Squint angle
FNBW = 1/16;        % First null beam width, only affects azimuth signal envelope
%% Target generatoin
r_target = [-25, -25, 25];                  % Range coordinates of targets
a_target = [-50, -50 * tan(Theta_sq), 0];      % Azimuth coordinates of targets

d0_target = r_target * tan(Theta_sq) - a_target;   % Azimuth distance from targets to beam center at t=0

figure; scatter(r_target, a_target, 'filled');
x = min(r_target):0; y = x * tan(Theta_sq);
title('Target locations and beam center at t=0');
hold on;
scatter(0, 0, 'r', 'filled');
plot(x, y, 'Color', 'red', 'LineWidth', 1);
clear('x', 'y');
%% Reflected radar signal generation
ta = 1/Fa * ((0:Na-1) - Na/2);              % Azimuth sampling times
R0 = Rc * cos(Theta_sq);                       % Minimum range of image center

% Range of targets (row: azimuth time, column: target)
R = sqrt((Rc * sin(Theta_sq) + a_target - V * ta').^2 + (R0 + r_target).^2);
d = d0_target + V * ta';

figure; plot(R, (0:Na-1) - Na/2);
title('Target range');

tr = 1/Fr * ((0:Nr-1) - Nr/2);              % Range sampling times
tr_real = tr + 2*Rc/c;                      % Actual time delay at sampling times
R = reshape(R, 256, 1, numel(r_target));    % Use dim3 to indicate targets, row: azimuth time
t_s = tr_real - R .* 2 ./ c;                % Time of reflected signals, row: azimuth time，column: range time

p_delay = exp(-1i * 4*pi * f0 * R ./ c);    % Phase produced by reflected signal delays
p_chirp = exp(1i * pi * Kr * (t_s).^2);     % Phase of LFM signal
wr = double(abs(t_s) <= Tr/2);              % Range Signal envelope
wa = sinc(2/FNBW * (reshape(d, 256, 1, numel(r_target)) * cos(Theta_sq) ./ R)).^2;   % Azimuth signal envelope
s = sum(wa .* wr .* p_chirp .* p_delay, 3); % Received signal, summed to all targets

clear('d', 'tr_real', 'R', 't_s', 'p_delay', 'p_chirp', 'wr', 'wa');
figure;
subplot(2, 2, 1); surf(real(s), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('(a)Real');
subplot(2, 2, 2); surf(imag(s), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('(b)Imag');
subplot(2, 2, 3); surf(abs(s), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('(c)Amplitude');
subplot(2, 2, 4); surf(angle(s), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('(d)Phase');
%% RDA
fr = Fr/Nr * ((0:Nr-1) - Nr/2);             % Range frequencies
% Range matched filter
Hr = exp(1i * pi / Kr * fr.^2);

[~, i] = max(abs(fr) <= Kr * Tr /2);
Interval = [i, i - 1 + sum(abs(fr) <= Kr * Tr /2)];
Hr = ifftshift(AddWindow(Hr, Interval));    % Kaiser window with length equal to bandwidth
clear('i', 'Interval');

S = fft(s, Nr, 2);
s_comp = ifft(Hr .* S, Nr, 2);              % Range matched filtering

clear('Hr', 'wk', 'S');

fa = Fa/Na * ((0:Na-1) - Na/2);             % Azimuth frequencies
fc = 2 * f0 * V * sin(Theta_sq) / c;        % Doppler center frequency
fa = round((fc - fa) / Fa) * Fa + fa;       % Shift azimuth frequency from [-Fa/2, Fa/2] to [fc-Fa/2, fc+Fa/2]
R0_r = (Rc + c * tr / 2) * cos(Theta_sq);   % Minimum range of range cells
D = sqrt(1 - (c * fa / 2 / f0 / V).^2);
Ha = exp(1i * 4*pi * f0 / c * R0_r .* D.'); % Accurate azimuth matched filter, column: R0
Ha = ifftshift(Ha, 1);                      % Shift zero frequency to left end

S_comp = fft(s_comp, Na);
clear('s_comp');

tc = 2/c * (R0_r ./ D.' - R0_r);            % Accurate Range cell migration delay
% tc = R0_r * c .* (fa').^2 / (4 * (f0 * V)^2);     % Approximate Range cell migration delay
tc = tc - 2/c * (Rc - R0);                  % Subtract time delay induced by squint angle

S_comp = fftshift(S_comp, 1);               % Shift to actual frequency
S_rcmc = zeros(Na, Nr);                     % Used to store calibrated signal

% Frequency modulation rate of secondary range compression, row: azimuth frequency
K_src = 2 * V^2 * f0^3 * (D.').^3 ./ (c * R0 .* (fa.').^2);

N_intp = 16;     % Number of points used in an interpolation
P_intp = 8;     % Precision of interpolation
n_intp = -N_intp/2 + 1:N_intp/2;            % Adjacent points to be used in interpolation
dt_intp = (0:P_intp-1) / P_intp;            % Discrete interpolation times
IK = sinc(dt_intp - n_intp');

n_c = Fr * (tr + tc);
n_c = round(n_c * P_intp) / P_intp;         % Interpolation time discretization
for i = 1:Na
    % SRC filter
    H_src = exp(-1i * pi * fr.^2 ./ K_src(i));
    H_src = ifftshift(H_src);
    S_src = ifft(fft(S_comp(i, :)) .* H_src);
    for j = 1:Nr
        S_local = zeros(1, N_intp);
        i_n = floor(n_c(i, j));
        d_n = n_c(i, j) - i_n;
        index = n_intp + i_n + Nr/2 + 1;   % Original data points used to interpolate the point (i,j)
        v = index >= 1 & index <= Nr;      % Nonzero elements in these data points
        S_local(v ~= 0) = S_src(index(v));
        S_rcmc(i, j) = S_local * IK(:, d_n * P_intp + 1);
    end
end

clear('S_comp', 'K_src', 'H_src', 's_src', 'S_src', 'S_cf', 'n_c', 'tc');
clear('n_intp', 'dt_intp', 'IK', 'S_local', 'i_n', 'd_n', 'index', 'v');
clear('i', 'j');

S_rcmc = ifftshift(S_rcmc, 1);

s_i = ifft(Ha .* S_rcmc, Na);              % Azimuth matched filtering

clear('Ha', 'S_rcmc');

figure;
surf(abs(s_i), 'EdgeColor', 'none');
view(2); xlim([0, 255]); ylim([0, 255]); title('Image (RDA)');
%% CSA
fa = Fa/Na * ((0:Na-1) - Na/2);
fc = 2 * f0 * V * sin(Theta_sq) / c;
fa = round((fc - fa) / Fa) * Fa + fa;
D = sqrt(1 - (c * fa / 2 / f0 / V).^2)';
D_ref = cos(Theta_sq);
Km = Kr ./ (1 - Kr * (c * R0 * (fa').^2)/(2 * V^2 * f0^3 * D_ref^3));

tr_rcmr = tr + 2 * Rc / c - 2 * R0 / c ./ D;    % Time in scaling signal
s_sc = exp(1i * pi * Km .* (D_ref./D - 1) .* tr_rcmr.^2);   % LFM signal used to perform chirp scaling
S_rcmr = fft(s) .* ifftshift(s_sc, 1);     % Residual RCMC
clear('tr_rcmr', 's_sc');

fr = Fr/Nr * ((0:Nr-1) - Nr/2);
S_rcmr = fft(S_rcmr, Nr, 2);
% Range compression
H_rcomp = exp(1i * pi * D ./ Km ./ D_ref .* fr.^2);
[~, i] = max(abs(fr) <= Kr * Tr /2);
Interval = [i, i - 1 + sum(abs(fr) <= Kr * Tr /2)];
H_rcomp = AddWindow(H_rcomp, Interval);
% Bulk RCMC
H_rcmb = exp(1i * 4*pi/c * R0 * (1./D - 1/D_ref) .* fr);
S_rcmc = ifft(S_rcmr .* ifftshift(H_rcomp .* H_rcmb), Nr, 2);

clear('S_rcmr', 'H_rcomp', 'i', 'Interval', 'H_rcmb')

R0_r = (Rc + c * tr / 2) * cos(Theta_sq);
Ha = exp(1i * 4*pi * f0 / c * R0_r .* D);   % Azimuth matched filter
% Remove additional phase introduced in CS process
p_r = exp(-1i * 4*pi/c^2 * Km .* (1 - D ./ D_ref) .* ((R0_r - R0) ./ D).^2);
s_i = ifft(S_rcmc .* ifftshift(Ha .* p_r, 1));

clear('Ha', 'p_r', 'S_rcmc');
figure;
surf(abs(s_i), 'EdgeColor', 'none');
view(2); xlim([0, 255]); ylim([0, 255]); title('Image (CSA)');
%% wkA
fa = Fa/Na * ((0:1:Na-1) - Na/2);
fr = Fr/Nr * ((0:1:Nr-1) - Nr/2);
fc = 2 * f0 * V * sin(Theta_sq) / c;
fa = round((fc - fa) / Fa) * Fa + fa;
% Reference function (calibration of a target locating at image center)
RF = (4*pi*R0/c) * sqrt((f0 + fr).^2 - c^2 * (fa').^2 / (4 * V^2)) + pi * fr.^2 / Kr;
RF = exp(1i * RF);
S = fftshift(fft2(s));
S_RFM = S .* RF;

clear('RF', 'S');

S_i = zeros(Na, Nr);
for i = 1:Na
    f_i = sqrt((f0 + fr).^2 - c^2 * fc^2 / 4 / V^2) - f0;
    f_i = (max(f_i) + min(f_i))/2 + fr;
    f_i = sqrt((f0 + f_i).^2 + c^2 * fa(i)^2 / 4 / V^2) - f0;
    f_i = f_i ./ (Fr/Nr);
    M = sinc(f_i - ((0:Nr-1) - Nr/2)');
    S_i(i, :) = S_RFM(i, :) * M;
end

s_i = ifft2(ifftshift(S_i .* exp(-1i * 4*pi*Rc/c * fr)));
clear('S_i', 'f_i', 'i', 'M', 'S_RFM');

figure;
surf(abs(s_i), 'EdgeColor', 'none');
view(2); xlim([0, 255]); ylim([0, 255]); title('Image (wkA)');