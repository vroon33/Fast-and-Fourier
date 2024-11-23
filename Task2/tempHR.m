% Load the E1.mat file
load('E1.mat'); % Ensure E1.mat contains a variable named 'E1'

% Sampling frequency
Fs = 128; % in Hz

% Calculate the time vector
t = (0:length(E1)-1) / Fs;

% Plot the ECG signal
figure;
plot(t, E1);
xlabel('Time (s)');
ylabel('ECG Amplitude');
title('Noise-Free ECG Signal');
grid on;

% Detect R-peaks in the ECG signal
[pks, locs] = findpeaks(E1, 'MinPeakHeight', 0.5, 'MinPeakDistance', 0.6 * Fs);

% Plot the R-peaks on the ECG signal
hold on;
plot(locs / Fs, pks, 'ro');
legend('ECG Signal', 'R-Peaks');

% Calculate time intervals between R-peaks (in seconds)
R_intervals = diff(locs) / Fs;

% Calculate BPM from the intervals
BPM = 60 ./ R_intervals;

% Average BPM
avg_BPM = mean(BPM);

% Display results
disp(['Estimated Average Heart Rate: ', num2str(avg_BPM), ' BPM']);
