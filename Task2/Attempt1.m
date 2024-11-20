% function analyzeECG()
%     % Load ECG signals from .mat files
%     E1_data = load('E1.mat');
%     E2_data = load('E2.mat');
%     E3_data = load('E3.mat');
    
%     % Extract the signal data
%     % Note: You might need to adjust the field names based on your .mat file structure
%     % Commonly the variable names might be 'data', 'signal', or the same as the filename
%     % Use whos(matfile('E1.mat')) to check the variable names in your .mat files
%     E1 = double(extractSignal(E1_data));
%     E2 = double(extractSignal(E2_data));
%     E3 = double(extractSignal(E3_data));
    
%     % Ensure signals are column vectors
%     E1 = E1(:);
%     E2 = E2(:);
%     E3 = E3(:);
    
%     % Sampling rate is given as 128 Hz
%     fs = 128;
    
%     % Process clean signal (E1)
%     hr_E1 = estimateHeartRate(E1, fs, false);
%     fprintf('E1 Heart Rates per minute: \n');
%     disp(hr_E1);
    
%     % Process noisy signals (E2 and E3) with noise removal
%     hr_E2 = estimateHeartRate(E2, fs, true);
%     fprintf('E2 Heart Rates per minute: \n');
%     disp(hr_E2);
    
%     hr_E3 = estimateHeartRate(E3, fs, true);
%     fprintf('E3 Heart Rates per minute: \n');
%     disp(hr_E3);
    
%     % Plot results
%     plotResults(hr_E1, hr_E2, hr_E3);
    
%     % Optionally visualize the signals
%     figure('Name', 'Raw ECG Signals');
%     subplot(3,1,1); plot(E1); title('E1 (Clean)'); xlabel('Samples'); ylabel('Amplitude');
%     subplot(3,1,2); plot(E2); title('E2 (Noisy)'); xlabel('Samples'); ylabel('Amplitude');
%     subplot(3,1,3); plot(E3); title('E3 (Noisy)'); xlabel('Samples'); ylabel('Amplitude');
% end

% function signal = extractSignal(matData)
%     % Get the fieldnames of the structure
%     fields = fieldnames(matData);
%     % Take the first field (assuming it contains the signal)
%     signal = matData.(fields{1});
% end

% function hr = estimateHeartRate(signal, fs, denoise)
%     if denoise
%         % Noise removal steps
%         signal = removeNoise(signal, fs);
%     end
    
%     % Calculate samples per minute
%     samples_per_minute = fs * 60;
    
%     % Initialize heart rate array
%     num_minutes = max(1, floor(length(signal) / samples_per_minute));
%     hr = zeros(num_minutes, 1);
    
%     for i = 1:num_minutes
%         % Get one minute of data
%         start_idx = (i-1) * samples_per_minute + 1;
%         end_idx = min(i * samples_per_minute, length(signal));
%         minute_data = signal(start_idx:end_idx);
        
%         % Detect R peaks and calculate heart rate
%         hr(i) = detectRPeaks(minute_data, fs);
%     end
% end

% function cleaned_signal = removeNoise(signal, fs)
%     % 1. Apply bandpass filter (5-15 Hz for QRS complex)
%     [b, a] = butter(4, [5 15]/(fs/2), 'bandpass');
%     filtered = filtfilt(b, a, signal);
    
%     % 2. Remove baseline wander using high-pass filter
%     [b, a] = butter(2, 0.5/(fs/2), 'high');
%     filtered = filtfilt(b, a, filtered);
    
%     % 3. Apply moving average to smooth the signal
%     window_size = round(0.05 * fs); % 50ms window
%     cleaned_signal = movmean(filtered, window_size);
% end

% function hr = detectRPeaks(signal, fs)
%     % Normalize signal
%     signal = signal - mean(signal);
%     if std(signal) ~= 0
%         signal = signal / std(signal);
%     end
    
%     % Find peaks
%     min_peak_distance = round(0.4 * fs); % Minimum 400ms between peaks
%     min_peak_height = 0.5; % Adjust based on signal characteristics
%     [pks, locs] = findpeaks(signal, 'MinPeakDistance', min_peak_distance, ...
%                          'MinPeakHeight', min_peak_height);
    
%     % Calculate heart rate (beats per minute)
%     if length(locs) > 1
%         % Average time between peaks
%         mean_interval = mean(diff(locs)) / fs;
%         hr = 60 / mean_interval;
%     else
%         hr = 0; % No peaks detected
%     end
% end

% function plotResults(hr_E1, hr_E2, hr_E3)
%     figure('Name', 'Heart Rate Analysis');
    
%     % Time vector (minutes)
%     t1 = 1:length(hr_E1);
%     t2 = 1:length(hr_E2);
%     t3 = 1:length(hr_E3);
    
%     % Plot heart rates
%     subplot(3,1,1);
%     plot(t1, hr_E1, 'b-o', 'LineWidth', 1.5);
%     title('Clean Signal (E1) Heart Rate');
%     xlabel('Time (minutes)');
%     ylabel('Heart Rate (BPM)');
%     grid on;
%     ylim([40 120]); % Typical heart rate range
    
%     subplot(3,1,2);
%     plot(t2, hr_E2, 'r-o', 'LineWidth', 1.5);
%     title('Denoised Signal (E2) Heart Rate');
%     xlabel('Time (minutes)');
%     ylabel('Heart Rate (BPM)');
%     grid on;
%     ylim([40 120]);
    
%     subplot(3,1,3);
%     plot(t3, hr_E3, 'g-o', 'LineWidth', 1.5);
%     title('Denoised Signal (E3) Heart Rate');
%     xlabel('Time (minutes)');
%     ylabel('Heart Rate (BPM)');
%     grid on;
%     ylim([40 120]);
% end

% function visualizeSignalProcessing(signal, fs)
%     % Create time vector
%     t = (0:length(signal)-1)/fs;
    
%     figure('Name', 'Signal Processing Steps');
    
%     % Original signal
%     subplot(4,1,1);
%     plot(t, signal);
%     title('Original Signal');
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     grid on;
    
%     % After bandpass filter
%     [b, a] = butter(4, [5 15]/(fs/2), 'bandpass');
%     filtered = filtfilt(b, a, signal);
%     subplot(4,1,2);
%     plot(t, filtered);
%     title('After Bandpass Filter (5-15 Hz)');
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     grid on;
    
%     % After baseline removal
%     [b, a] = butter(2, 0.5/(fs/2), 'high');
%     baseline_removed = filtfilt(b, a, filtered);
%     subplot(4,1,3);
%     plot(t, baseline_removed);
%     title('After Baseline Removal');
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     grid on;
    
%     % Final cleaned signal with detected peaks
%     window_size = round(0.05 * fs);
%     cleaned = movmean(baseline_removed, window_size);
%     subplot(4,1,4);
%     plot(t, cleaned);
%     hold on;
    
%     % Find and plot peaks
%     [pks, locs] = findpeaks(cleaned, 'MinPeakDistance', round(0.4 * fs), ...
%                            'MinPeakHeight', 0.5 * std(cleaned));
%     plot(t(locs), pks, 'ro', 'MarkerSize', 10);
%     title('Cleaned Signal with Detected R-Peaks');
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     grid on;
%     hold off;
% end

% RUn the FFT of E1, E2, and E3
% function analyzeFFT()
%     % Load .mat files
%     E1 = load('E1.mat');
%     E2 = load('E2.mat');
%     E3 = load('E3.mat');
    
%     % Get the first field name from each structure
%     E1_vars = fieldnames(E1);
%     E2_vars = fieldnames(E2);
%     E3_vars = fieldnames(E3);
    
%     % Get the actual data
%     signal1 = E1.(E1_vars{1});
%     signal2 = E2.(E2_vars{1});
%     signal3 = E3.(E3_vars{1});
    
%     % Get sampling frequency
%     fs = 128; % Update if needed
    
%     % Compute FFT for each signal
%     fft1 = fft(signal1);
%     fft2 = fft(signal2);
%     fft3 = fft(signal3);
    
%     % Calculate frequency axis
%     L = length(signal1);
%     f = fs * (0:(L/2))/L;
    
%     % Get single-sided spectrum
%     P1 = abs(fft1(1:L/2+1));
%     P2 = abs(fft2(1:L/2+1));
%     P3 = abs(fft3(1:L/2+1));
    
%     % Plot results
%     figure('Position', [100 100 800 600]);
    
%     subplot(3,1,1)
%     plot(f, P1)
%     title('FFT of E1')
%     xlabel('Frequency (Hz)')
%     ylabel('Magnitude')
%     grid on
    
%     subplot(3,1,2)
%     plot(f, P2)
%     title('FFT of E2')
%     xlabel('Frequency (Hz)')
%     ylabel('Magnitude')
%     grid on
    
%     subplot(3,1,3)
%     plot(f, P3)
%     title('FFT of E3')
%     xlabel('Frequency (Hz)')
%     ylabel('Magnitude')
%     grid on
    
%     sgtitle('FFT Analysis of Signals')
% end

% Load the .mat file and check its contents
% load('E1.mat');
% % vars = whos('-file', 'E1.mat');
% % disp('Variables in the E1.mat file:');
% % disp(vars);

% % Get the first variable from the .mat file (assuming it's the ECG signal)
% data_fields = fieldnames(load('E1.mat'));
% ecg_signal = getfield(load('E1.mat'), data_fields{1});

% % Get the sampling frequency
% Fs = 128;  % Given as 128 samples/second


% % Function to estimate heart rate per minute
% function hr = estimate_heart_rate(signal, window_samples, Fs)
%     % Find peaks (R peaks in ECG)
%     [peaks, locs] = findpeaks(signal, 'MinPeakHeight', mean(signal) + std(signal), ...
%                             'MinPeakDistance', 0.5*Fs);  % Minimum 0.5s between peaks
    
%     seconds_per_minute = 60;

%     % Calculate heart rate (beats per minute)
%     time_diff = diff(locs)/Fs;  % Time between peaks in seconds
%     inst_hr = seconds_per_minute./time_diff;  % Instantaneous heart rate
    
%     % Average heart rate for the window
%     hr = mean(inst_hr);
% end

% % Process signal minute by minute
% window_duration = 60;  % 1 minute
% window_samples = window_duration * Fs;
% total_minutes = floor(length(ecg_signal)/window_samples);

% % Initialize arrays for time and heart rate
% time = zeros(1, total_minutes);
% heart_rate = zeros(1, total_minutes);

% % Calculate heart rate for each minute
% for i = 1:total_minutes
%     start_idx = (i-1)*window_samples + 1;
%     end_idx = i*window_samples;
    
%     window_signal = ecg_signal(start_idx:end_idx);
%     heart_rate(i) = estimate_heart_rate(window_signal, window_samples, Fs);
%     time(i) = i;
% end

% % Plot results
% figure;
% plot(time, heart_rate, 'b-o', 'LineWidth', 2);
% title('Heart Rate Estimation Over Time');
% xlabel('Time (minutes)');
% ylabel('Heart Rate (BPM)');
% grid on;
