% Read the audio file
[y, Fs] = audioread('Bird1.wav');

% y = y(0.95*Fs : 2.5*Fs);
% Compute FFT
y = [y zeros(1, 2^nextpow2(length(y)) - length(y))];
N = length(y);
Y = fft(y);
Y_shifted = fftshift(Y);
mag_spectrum = abs(Y_shifted);  % Removed /N scaling

% Create frequency vector centered at 0
freq = (-N/2:N/2-1)*(Fs/N);

% Only consider the positive half of the spectrum
positive_freq_idx = freq > 0;

positive_mag_spectrum = mag_spectrum(positive_freq_idx);
positive_freq = freq(positive_freq_idx);

% Find peaks in the positive frequency spectrum
[peaks, locs] = findpeaks(positive_mag_spectrum, 'MinPeakProminence', max(positive_mag_spectrum)/20);
peak_freqs = positive_freq(locs);

% Sort peaks by magnitude to find dominant frequencies
[sorted_peaks, idx] = sort(peaks, 'descend');
dominant_freqs = peak_freqs(idx);

% Limit to the first 5 positive peaks
top_peaks = dominant_freqs(1:min(5, length(dominant_freqs)));

% Noise analysis
% Estimate noise floor (mean of bottom 10% of magnitudes)
sorted_mags = sort(mag_spectrum);
noise_floor = mean(sorted_mags(1:floor(length(sorted_mags)*0.1)));

% Calculate Signal-to-Noise Ratio (SNR)
signal_power = mean(abs(y).^2);
noise_power = noise_floor^2;
SNR = 10*log10(signal_power/noise_power);

% Plot frequency spectrum with peaks marked
figure;

% Main spectrum plot (only positive frequencies)
plot(positive_freq, positive_mag_spectrum);
hold on;
% Plot detected peaks
plot(top_peaks, sorted_peaks(1:length(top_peaks)), 'r*', 'MarkerSize', 10);
% Plot noise floor
plot(freq, noise_floor*ones(size(freq)), 'r--', 'LineWidth', 1);

title('Frequency Spectrum with Peak Detection (Positive Frequencies)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs / 2]);
grid on;
set(gca, 'YScale', 'log');
legend('Spectrum', 'Peak Frequencies', 'Noise Floor');

% Print analysis results
fprintf('\nSignal Analysis Results:\n');
fprintf('------------------------\n');
fprintf('Sampling Frequency: %d Hz\n', Fs);
fprintf('Signal Duration: %.2f seconds\n', N/Fs);
fprintf('Estimated SNR: %.2f dB\n', SNR);
fprintf('Noise Floor Level: %.2e\n', noise_floor);

fprintf('\nDominant Frequencies (top 5):\n');
for i = 1:min(5,length(dominant_freqs))
    fprintf('Peak %d: %.1f Hz (Magnitude: %.2e)\n', i, dominant_freqs(i), sorted_peaks(i));
end

% Plot the signal waveform
figure(2);
t = (0:N-1)/Fs;
plot(t, y);
title('Signal Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;