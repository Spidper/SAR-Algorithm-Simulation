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
Theta = 22.5 * pi/180;  % Squint angle
FNBW = 1/16;        % First null beam width, only affects azimuth signal envelope
%% Target generatoin
r_target = [-25, -25, 25];                  % Range coordinates of targets
a_target = [-50, -50 * tan(Theta), 0];      % Azimuth coordinates of targets

d0_target = r_target * tan(Theta) - a_target;   % Azimuth distance from targets to beam center at t=0

figure; scatter(r_target, a_target, 'filled');
xlim([-50, 50]); ylim([-75, 25]);
x = -50:0; y = x * tan(Theta);
title('Target locations and beam center at t=0');
hold on;
scatter(0, 0, 'r', 'filled');
plot(x, y, 'Color', 'red', 'LineWidth', 1);
%% Reflected radar signal generation
ta = 1/Fa * ((0:Na-1) - Na/2);              % Azimuth sampling times
R0 = Rc * cos(Theta);                       % Minimum range of image center

% Range of targets (row: azimuth time, column: target)
R = sqrt((Rc * sin(Theta) + a_target - V * ta').^2 + (R0 + r_target).^2);
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
wa = sinc(2/FNBW * (reshape(d, 256, 1, numel(r_target)) * cos(Theta) ./ R)).^2;   % Azimuth signal envelope
s = sum(wa .* wr .* p_chirp .* p_delay, 3); % Received signal, summed to all targets

figure;
subplot(2, 2, 1); surf(real(s), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('(a)Real');
subplot(2, 2, 2); surf(imag(s), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('(b)Imag');
subplot(2, 2, 3); surf(abs(s), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('(c)Amplitude');
subplot(2, 2, 4); surf(angle(s), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('(d)Phase');
%% Pulse compression
fr = Fr/Nr * ((0:Nr-1) - Nr/2);             % Range frequencies
% Range matched filter
Hr = exp(1i * pi / Kr * fr.^2);
Hr = ifftshift(Hr .* (abs(fr) <= Kr * Tr /2));  % Shift zero frequency to left end

wk = kaiser(sum(abs(fr) <= Kr * Tr /2), 2.5);
[~, i] = max(abs(fr) <= Kr * Tr /2);
wk = [zeros(1, i-1), wk'];
wk = [wk, zeros(1, Nr - numel(wk))];        % Kaiser window with length equal to bandwidth
Hr = Hr .* ifftshift(wk);

S = fft(s, Nr, 2);
s_c = ifft(Hr .* S, Nr, 2);                 % Range matched filtering

figure;
surf(abs(s_c), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('Signal amplitude after pulse compression');
%% Direct azimuth matched filtering
fa = Fa/Na * ((0:Na-1) - Na/2);             % Azimuth frequencies
fc = 2 * f0 * V * sin(Theta) / c;           % Doppler center frequency
fa = round((fc - fa) / Fa) * Fa + fa;       % Shift azimuth frequency from [-Fa/2, Fa/2] to [fc-Fa/2, fc+Fa/2]
R0_r = (Rc + c * tr / 2) * cos(Theta);      % Minimum range of sampling points
D = sqrt(1 - (c * fa / 2 / f0 / V).^2);
% Ka = 2 * V.^2 * f0 / c / R0;              % Azimuth frequency modulation rate
% Ha = exp(-1i * pi / Ka * fa.^2).';        % Approximate azimuth matched filter
Ha = exp(1i * 4*pi * f0 / c * R0_r .* D.'); % Accurate azimuth matched filter, column: R0
Ha = ifftshift(Ha, 1);                      % Shift zero frequency to left end

S_c = fft(s_c, Na);
surf(abs(S_c), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('Azimuth spectrum amplitude');

s_i = ifft(Ha .* S_c, Na);                  % Azimuth matched filtering without RCMC
figure;
surf(abs(s_i), 'EdgeColor', 'none');
view(2); xlim([0, 255]); ylim([0, 255]); title('Signal amplitude after azimuth matched filtering');
%% Secondary range compression, range cell migration calibration and azimuth matched filtering
tc = 2/c * (R0_r ./ D.' - R0_r);            % Accurate Range cell migration delay
% tc = R0_r * c .* (fa').^2 / (4 * (f0 * V)^2);     % Approximate Range cell migration delay
tc = tc - 2/c * (Rc - R0);                  % Subtract time delay induced by squint angle

% Frequency modulation rate of secondary range compression, row: azimuth frequency, column: R0
K_src = 2 * V^2 * f0^3 * (D.').^3 ./ (c * R0_r .* (fa.').^2);

S_c = fftshift(S_c, 1);                     % Shift to actual frequency
S_cc = zeros(Na, Nr);                       % Used to store calibrated signal
for i = 1:Na
    n_c = Fr * (tr + tc(i, :));             % Interpolate with sinc(Fr*t), y(t)=sum(x(n)sinc(Fr*t-n))
    M = sinc(n_c - ((0:Nr-1) - Nr/2)');     % Interpolation matrix
    
    % SRC filter, column: R0, because required filter changes with R0
    H_src = exp(-1i * pi * (fr').^2 ./ K_src(i, :));
    H_src = ifftshift(H_src, 1);            % Shift zero frequency to left end
    
    S_cf = fft(S_c(i, :));                  % Range spectrum at i-th azimuth frequency
    % Filter the range spectrum using filter corresponding to j-th sampling time,
    % use diag() to extract filtered signal at j-th sampling time
    S_src = diag(ifft((S_cf.') .* H_src)).';
    
    S_cc(i, :) = S_src * M;                 % Store filtered signal into i-th row
end
S_c = ifftshift(S_c, 1);

surf(abs(S_cc), 'EdgeColor', 'none');
view(2); xlim([0,255]); ylim([0,255]); title('Spectrum amplitude after SRC and RCMC');

S_cc = ifftshift(S_cc, 1);
s_ic = ifft(Ha .* S_cc, Na);                 % Azimuth matched filtering after RCMC
figure;
surf(abs(s_ic), 'EdgeColor', 'none');
view(2); xlim([0, 255]); ylim([0, 255]); title('Signal amplitude after azimuth filtering with SRC and RCMC');