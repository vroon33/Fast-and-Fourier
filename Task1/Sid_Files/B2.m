% Read the audio file
[y, Fs] = audioread('./Fast-and-Fourier/Task1/bird2.wav');

y = noisefilter(y, Fs);
sound(y, Fs);

% Compute FFT
N = length(y);
fft_length = 2 * Fs;  % Match the FFT length
Y = fft(y, fft_length)';
Y_shifted = fftshift(Y);  % Center the FFT
mag_spectrum = abs(Y_shifted);  % Get magnitude spectrum

% Create frequency vector centered at 0
freq = (-fft_length/2:fft_length/2-1) * (Fs/fft_length);

% Find peaks with a minimum distance of 100 Hz
min_distance = 100 / (Fs / length(freq));  % Convert 100 Hz to index distance
[peaks, peak_locs] = findpeaks(mag_spectrum, 'MinPeakDistance', min_distance);

% Sort peaks by magnitude (descending order) and get the top 5
[sorted_peaks, sorted_indices] = sort(peaks, 'descend');
top_peaks = sorted_peaks(1:min(20, length(sorted_peaks)));  % Take top 10 peaks
top_peak_locs = peak_locs(sorted_indices(1:min(20, length(sorted_peaks))));
top_peak_freqs = freq(top_peak_locs);  % Corresponding frequencies

% Plot the spectrum and highlight the top 5 peaks
figure;
plot(freq, mag_spectrum);
hold on;
plot(top_peak_freqs, top_peaks, 'r*', 'MarkerSize', 8);  % Highlight top peaks
title('Frequency Spectrum with Top 10 Peaks');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs / 2]);
grid on;
% set(gca, 'YScale', 'log');

% Display the top 10 peaks and their frequencies
fprintf('Top 10 Peaks:\n');
fprintf('Frequency (Hz): %.2f\tMagnitude: %.2f\n', [top_peak_freqs(2:2:end); top_peaks(2:2:end)]);

% Print audio info
fprintf('Sampling Frequency: %d Hz\n', Fs);
fprintf('Signal Duration: %.2f seconds\n', N/Fs);