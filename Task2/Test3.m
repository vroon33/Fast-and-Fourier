clc, clearvars, close all;

% Load data and set sampling frequency
load('E3.mat');
Fs = 128;

% Apply noise filtering with multiple stages
% First notch filter at 22 Hz with wider bandwidth
wo1 = 22/(Fs/2);
bw1 = wo1/100;  % Much wider bandwidth for very gentle filtering
[b1, a1] = iirnotch(wo1, bw1);
E3 = filtfilt(b1, a1, E3);

% Second notch filter at 50 Hz with wider bandwidth
wo2 = 50/(Fs/2);
bw2 = wo2/100;  % Much wider bandwidth for very gentle filtering
[b2, a2] = iirnotch(wo2, bw2);
E3 = filtfilt(b2, a2, E3);

% Create time vector
t = (0:length(E3)-1)/Fs;

% Adaptive threshold for R-peak detection
window_duration = 5;
window_samples = window_duration * Fs;
threshold = zeros(size(E3));

for i = 1:window_samples:length(E3)
    window_end = min(i + window_samples - 1, length(E3));
    window_data = E3(i:window_end);
    threshold(i:window_end) = mean(window_data) + 1.5*std(window_data);
end

% R-peak detection
[~, R_locs] = findpeaks(E3, ...
    'MinPeakHeight', mean(threshold), ...
    'MinPeakDistance', 0.25*Fs);

% Estimate BPM for each second
bpmVals = zeros(1, floor(length(E3)/Fs));

for sec = 1:length(bpmVals)
    % Define time window (30 seconds before and after)
    start_time = max(1, (sec-1)*Fs - 30*Fs);
    end_time = min(length(E3), (sec-1)*Fs + 30*Fs);
    
    % Find R-peaks in the window
    window_peaks = R_locs((R_locs >= start_time) & (R_locs <= end_time));
    
    % Calculate BPM
    if length(window_peaks) > 1
        RR_intervals = diff(window_peaks) / Fs;
        bpmVals(sec) = 60 / mean(RR_intervals);
    else
        if sec > 1
            bpmVals(sec) = bpmVals(sec-1);
        else
            bpmVals(sec) = 0;
        end
    end
end

% Visualization
figure;
subplot(2,1,1);
plot(t, E3);
hold on;
plot(t(R_locs), E3(R_locs), 'ro', 'MarkerSize', 8);
xlabel('Time (s)');
ylabel('ECG Amplitude');
title('ECG Signal with R-peaks');
grid on;

subplot(2,1,2);
plot(t(1:length(bpmVals)), bpmVals, 'b-');
xlabel('Time (s)');
ylabel('Heart Rate (BPM)');
title('Heart Rate Over Time');
grid on;

% Summary statistics
avgBPM = mean(bpmVals);
fprintf('Average heart rate: %.1f BPM\n', avgBPM);
fprintf('Min heart rate: %.1f BPM\n', min(bpmVals));
fprintf('Max heart rate: %.1f BPM\n', max(bpmVals));

% ... (keep all your previous code until the summary statistics) ...

% Calculate how many complete minutes of data we want
samples_per_minute = 60 * Fs;
total_complete_minutes = floor(length(E3)/(samples_per_minute));
samples_to_keep = total_complete_minutes * samples_per_minute;

% Trim the signal to exact multiple of minutes
E3_trimmed = E3(1:samples_to_keep);

% Create time vector for trimmed signal
t = (0:length(E3_trimmed)-1)/Fs;

% Figure 1: Plot full trimmed ECG signal
figure;
subplot(2,1,1);
plot(t, E3_trimmed);
xlabel('Time (s)');
ylabel('ECG Amplitude');
title('Noise-Free ECG Signal (Trimmed to Complete Minutes)');
grid on;

% Improved R-peak detection
% First find all R-peaks at once with adaptive threshold
window_duration = 5; % 5 seconds for calculating local statistics
window_samples = window_duration * Fs;
signal_length = length(E3_trimmed);
threshold = zeros(size(E3_trimmed));

% Calculate adaptive threshold
for i = 1:window_samples:signal_length
    window_end = min(i + window_samples - 1, signal_length);
    window_data = E3_trimmed(i:window_end);
    threshold(i:window_end) = mean(window_data) + 1.5*std(window_data);
end

% Find R-peaks using the adaptive threshold
[pks, R_locs] = findpeaks(E3_trimmed, ...
    'MinPeakHeight', mean(threshold), ...
    'MinPeakDistance', 0.25*Fs);  % Minimum 0.5 seconds between peaks

% Calculate BPM values for each minute
total_minutes = total_complete_minutes;
bpmVals = zeros(1, total_minutes);

% Overlay R-peak detection on the ECG plot
subplot(2,1,1);
hold on;
plot(t(R_locs), E3_trimmed(R_locs), 'ro', 'MarkerSize', 8);
title(sprintf('ECG Signal with Detected R-peaks (Total: %d)', length(R_locs)));

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

% Figure 1 (continued): Plot heart rate
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

% Add new code for plotting individual minutes
minutes_to_plot = total_minutes;
samples_per_minute = 60 * Fs;

for minute = 1:minutes_to_plot
    % Create a new figure for each minute
    figure;
    
    % Calculate time range for this minute
    start_sample = (minute-1) * samples_per_minute + 1;
    end_sample = minute * samples_per_minute;
    time_range = t(start_sample:end_sample);
    
    % Plot ECG signal for this minute
    plot(time_range, E3_trimmed(start_sample:end_sample), 'b');
    hold on;
    
    % Find and plot R-peaks for this minute
    minute_peaks = R_locs((R_locs >= start_sample) & (R_locs <= end_sample));
    plot(t(minute_peaks), E3_trimmed(minute_peaks), 'ro', 'MarkerSize', 8);
    
    % Add labels and title
    xlabel('Time (s)');
    ylabel('ECG Amplitude');
    title(sprintf('Minute %d - ECG Signal with R-peaks (BPM: %.1f)', ...
        minute, bpmVals(minute)));
    grid on;
    
    % Add text box with statistics for this minute
    stats_text = sprintf('R-peaks detected: %d\nBPM: %.1f', ...
        length(minute_peaks), bpmVals(minute));
    annotation('textbox', [0.7 0.7 0.2 0.2], ...
        'String', stats_text, ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'white');
    
    % Set consistent y-axis limits across all minutes
    ylim([min(E3_trimmed) max(E3_trimmed)]);
    
    % Adjust figure size and position
    set(gcf, 'Position', [100 100 800 400]);
end