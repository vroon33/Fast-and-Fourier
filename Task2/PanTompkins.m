load('E1.mat');
Fs = 128;

% Calculate how many complete minutes of data we want
samples_per_minute = 60*Fs;
total_complete_minutes = floor(length(E1)/(samples_per_minute));
samples_to_keep = total_complete_minutes * samples_per_minute;

% Trim the signal to exact multiple of minutes
ecg = E1(1:samples_to_keep);

% Create time vector for trimmed signal
t = (0:length(ecg)-1)/Fs;

% Step 1: Bandpass Filter (5-15 Hz)
[b, a] = butter(2, [5 15]/(Fs/2), 'bandpass');
ecg_filtered = filtfilt(b, a, ecg);

% Step 2: Derivative
derivative = zeros(size(ecg_filtered));
derivative(1:end-1) = diff(ecg_filtered);
derivative(end) = derivative(end-1);

% Step 3: Squaring
squared = derivative.^2;

% Step 4: Moving Window Integration
window_size = round(0.15*Fs); % 150ms window
moving_avg = movmean(squared, window_size);

% Adaptive thresholding
window_size_adapt = 2*Fs; % 2 second window for adaptive threshold
threshold = zeros(size(moving_avg));
for i = 1:window_size_adapt:length(moving_avg)
    window_end = min(i + window_size_adapt - 1, length(moving_avg));
    window_data = moving_avg(i:window_end);
    threshold(i:window_end) = mean(window_data) + 1.5*std(window_data);
end

% Find peaks with more sensitive parameters
min_distance = round(0.4*Fs); % Minimum 400ms between peaks (assuming max 150 BPM)
[pks, locs] = findpeaks(moving_avg, 'MinPeakHeight', mean(threshold), 'MinPeakDistance', min_distance);

% Find corresponding R-peaks in original ECG
R_peaks = zeros(size(locs));
window = round(0.2*Fs); % 100ms window to search for R peak

for i = 1:length(locs)
    % Define search window
    start_idx = max(1, locs(i) - window);
    end_idx = min(length(ecg), locs(i) + window);
    
    % Find maximum in original ECG within window
    [~, max_idx] = max(ecg(start_idx:end_idx));
    R_peaks(i) = start_idx + max_idx - 1;
end

% Calculate BPM for each minute
bpm = zeros(1, total_complete_minutes);
for minute = 1:total_complete_minutes
    start_sample = (minute-1) * samples_per_minute + 1;
    end_sample = minute * samples_per_minute;
    
    % Find R-peaks within this minute
    peaks_in_minute = sum(R_peaks >= start_sample & R_peaks <= end_sample);
    bpm(minute) = peaks_in_minute;
end

% Plot results
figure;
subplot(4,1,1);
plot(t, ecg);
hold on;
plot(t(R_peaks), ecg(R_peaks), 'ro');
title('Original ECG with R-peaks');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,2);
plot(t, ecg_filtered);
title('Bandpass Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,3);
plot(t, moving_avg);
title('Moving Average with Adaptive Threshold');
hold on;
plot(t, threshold, 'r--');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,4);
plot(1:total_complete_minutes, bpm, 'bo-');
title('Heart Rate by Minute');
xlabel('Minute');
ylabel('BPM');

% Print statistics
fprintf('\nECG Analysis Summary:\n');
fprintf('Total duration: %d complete minutes\n', total_complete_minutes);
fprintf('Total R-peaks detected: %d\n', length(R_peaks));
fprintf('Average R-peaks per minute: %.1f\n', mean(bpm));
fprintf('Min heart rate: %.1f BPM\n', min(bpm));
fprintf('Max heart rate: %.1f BPM\n', max(bpm));

% Plot individual minutes
for minute = 1:total_complete_minutes
    figure;
    
    % Calculate time range for this minute
    start_sample = (minute-1) * samples_per_minute + 1;
    end_sample = minute * samples_per_minute;
    time_range = t(start_sample:end_sample);
    
    % Plot ECG signal for this minute
    plot(time_range, ecg(start_sample:end_sample));
    hold on;
    
    % Find and plot R-peaks for this minute
    minute_peaks = R_peaks((R_peaks >= start_sample) & (R_peaks <= end_sample));
    plot(t(minute_peaks), ecg(minute_peaks), 'ro', 'MarkerSize', 8);
    
    % Add labels and title
    xlabel('Time (s)');
    ylabel('ECG Amplitude');
    title(sprintf('Minute %d - ECG Signal with R-peaks (BPM: %.1f)', ...
        minute, bpm(minute)));
    grid on;
    
    % Add text box with statistics
    stats_text = sprintf('R-peaks detected: %d\nBPM: %.1f', ...
        length(minute_peaks), bpm(minute));
    annotation('textbox', [0.7 0.7 0.2 0.2], ...
        'String', stats_text, ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'white');
end