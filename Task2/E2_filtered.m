% Load ECG signal
load('E2.mat');  % Load your noisy ECG file
load('E1.mat');
Fs = 128;

% Calculate how many complete minutes of data we want
samples_per_minute = 60 * Fs;
total_complete_minutes = floor(length(E2)/(samples_per_minute));
samples_to_keep = total_complete_minutes * samples_per_minute;

% Trim the signal to exact multiple of minutes
E2_trimmed = E2(1:samples_to_keep);

% Create time vector for trimmed signal
t = (0:length(E2_trimmed)-1)/Fs;

% Apply noise filtering with adjusted parameters
[b, a] = butter(3, 30/(Fs/2), 'low');  % Slightly reduced cutoff frequency
ecg_filtered = filtfilt(b, a, E2_trimmed);

% Normalize filtered signal to match E1's amplitude range
E1_trimmed = E1(1:samples_to_keep);  % Trim E1 to same length
E1_range = max(E1_trimmed) - min(E1_trimmed);
filtered_range = max(ecg_filtered) - min(ecg_filtered);
scaling_factor = E1_range / filtered_range;

% Scale the filtered signal with slightly adjusted scaling
scaling_factor = scaling_factor * 0.95;  % Reduce scaling slightly
ecg_filtered = ecg_filtered * scaling_factor;

% Adjust DC offset to match E1's mean
ecg_filtered = ecg_filtered - mean(ecg_filtered) + mean(E1_trimmed);

% Plot comparison of signals
figure;
subplot(3,1,1);
plot(t, E1_trimmed);
title('Reference Clean ECG (E1)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, E2_trimmed);
title('Original Noisy ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t, ecg_filtered);
title('Filtered ECG Signal (Amplitude Matched)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Normalize E1 and E2 signals
E1_normalized = E1 / max(abs(E1)); % Normalize E1 to range [-1, 1]
E2_normalized = E2 / max(abs(E2)); % Normalize E2 to range [-1, 1]

% Compute FFT of normalized E1
L1 = length(E1_normalized);
E1_fft = fft(E1_normalized);
P2_E1 = abs(E1_fft / L1); 
P1_E1 = P2_E1(1:L1/2+1);
P1_E1(2:end-1) = 2 * P1_E1(2:end-1);
f1 = Fs * (0:(L1/2)) / L1;

% Compute FFT of filtered E2
L2 = length(ecg_filtered);
E2_fft = fft(ecg_filtered);
P2_E2 = abs(E2_fft / L2);
P1_E2 = P2_E2(1:L2/2+1);
P1_E2(2:end-1) = 2 * P1_E2(2:end-1);
f2 = Fs * (0:(L2/2)) / L2;

% Plot FFTs
figure;
subplot(2,1,1);
plot(f1, P1_E1, 'b');
title('FFT of Normalized E1');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(f2, P1_E2, 'r');
title('FFT of Filtered E2');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Adjust figure layout
set(gcf, 'Position', [100, 100, 800, 600]);

% Improved R-peak detection on filtered signal
window_duration = 5; % 5 seconds for calculating local statistics
window_samples = window_duration * Fs;
signal_length = length(ecg_filtered);
threshold = zeros(size(ecg_filtered));

% Calculate adaptive threshold
for i = 1:window_samples:signal_length
    window_end = min(i + window_samples - 1, signal_length);
    window_data = ecg_filtered(i:window_end);
    threshold(i:window_end) = mean(window_data) + 1.5*std(window_data);
end

% Find R-peaks using the adaptive threshold
[pks, R_locs] = findpeaks(ecg_filtered, ...
    'MinPeakHeight', mean(threshold), ...
    'MinPeakDistance', 0.25*Fs);  % Minimum 0.25 seconds between peaks

% Calculate BPM values for each minute
total_minutes = total_complete_minutes;
bpmVals = zeros(1, total_minutes);

% Create new figure for R-peak detection and heart rate
figure;
subplot(2,1,1);
plot(t, ecg_filtered);
hold on;
plot(t(R_locs), ecg_filtered(R_locs), 'ro', 'MarkerSize', 8);
title(sprintf('Filtered ECG Signal with Detected R-peaks (Total: %d)', length(R_locs)));
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

for minute = 1:total_minutes
    % Time window for current minute
    start_sample = (minute-1) * 60 * Fs + 1;
    end_sample = minute * 60 * Fs;
    
    % Find R-peaks within this minute
    minute_peaks = R_locs((R_locs >= start_sample) & (R_locs <= end_sample));
    
    % Calculate BPM for this minute
    if length(minute_peaks) > 1
        RR_intervals = diff(minute_peaks) / Fs;  % Convert to seconds
        avg_RR = mean(RR_intervals);
        bpmVals(minute) = 60 / avg_RR;
        
        % Print number of R-peaks detected in this minute
        fprintf('Minute %d: %d R-peaks detected, BPM = %.1f\n', ...
            minute, length(minute_peaks), bpmVals(minute));
    else
        if minute > 1
            bpmVals(minute) = bpmVals(minute-1);
            fprintf('Warning: Less than 2 R-peaks detected in minute %d, using previous value\n', minute);
        else
            bpmVals(minute) = 0;
            fprintf('Warning: Less than 2 R-peaks detected in minute %d\n', minute);
        end
    end
end

% Calculate average BPM
avgBPM = mean(bpmVals);

% Plot heart rate
subplot(2,1,2);
time_mins = 1:total_minutes;
plot(time_mins, bpmVals, 'b-o', 'LineWidth', 2);
hold on;
plot([1 total_minutes], [avgBPM avgBPM], 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time (minutes)');
ylabel('Heart Rate (BPM)');
title('Heart Rate Over Time');
legend('Heart Rate', sprintf('Average: %.1f BPM', avgBPM));

% Adjust figure layout
set(gcf, 'Position', [100 100 800 600]);
subplot(2,1,1); box on;
subplot(2,1,2); box on;

% Print summary statistics
fprintf('\nECG Analysis Summary:\n');
fprintf('Total duration: %d complete minutes\n', total_minutes);
fprintf('Total R-peaks detected: %d\n', length(R_locs));
fprintf('Average R-peaks per minute: %.1f\n', length(R_locs)/total_minutes);
fprintf('Average heart rate: %.1f BPM\n', avgBPM);
fprintf('Min heart rate: %.1f BPM\n', min(bpmVals));
fprintf('Max heart rate: %.1f BPM\n', max(bpmVals));