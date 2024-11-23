% Load and plot the ECG signals
load('E1.mat'); % Ensure E1.mat contains variable E1
load('E2.mat'); % Ensure E2.mat contains variable E2
load('E3.mat'); % Ensure E3.mat contains variable E3

Fs = 128;       % Sampling frequency
T = 1 / Fs;     % Sampling period
L1 = length(E1); % Length of E1
L2 = length(E2); % Length of E2
L3 = length(E3); % Length of E3

t1 = (0:L1-1) * T; % Time vector for E1
t2 = (0:L2-1) * T; % Time vector for E2
t3 = (0:L3-1) * T; % Time vector for E3

% Compute FFT for E1
E1_fft = fft(E1);
P2_E1 = abs(E1_fft / L1); 
P1_E1 = P2_E1(1:L1/2+1);
P1_E1(2:end-1) = 2 * P1_E1(2:end-1);
f1 = Fs * (0:(L1/2)) / L1;

% Compute FFT for E2
E2_fft = fft(E2);
P2_E2 = abs(E2_fft / L2); 
P1_E2 = P2_E2(1:L2/2+1);
P1_E2(2:end-1) = 2 * P1_E2(2:end-1);
f2 = Fs * (0:(L2/2)) / L2;

% Compute FFT for E3
E3_fft = fft(E3);
P2_E3 = abs(E3_fft / L3); 
P1_E3 = P2_E3(1:L3/2+1);
P1_E3(2:end-1) = 2 * P1_E3(2:end-1);
f3 = Fs * (0:(L3/2)) / L3;

% Plot signals and their FFTs
figure;

% ECG Signals
subplot(3, 2, 1);
plot(t1, E1);
title('ECG Signal E1');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(3, 2, 3);
plot(t2, E2);
title('ECG Signal E2');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(3, 2, 5);
plot(t3, E3);
title('ECG Signal E3');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% FFT Spectrums
subplot(3, 2, 2);
plot(f1, P1_E1);
title('FFT of E1');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
grid on;

subplot(3, 2, 4);
plot(f2, P1_E2);
title('FFT of E2');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
grid on;

subplot(3, 2, 6);
plot(f3, P1_E3);
title('FFT of E3');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
grid on;
