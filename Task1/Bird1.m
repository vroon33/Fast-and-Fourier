% Read the audio file
[y, Fs] = audioread('bird1.wav');

% Compute FFT
N = length(y);
Y = fft(y);
Y_shifted = fftshift(Y);  % Center the FFT
mag_spectrum = abs(Y_shifted);  % Get magnitude spectrum

% Create frequency vector centered at 0
freq = (-N/2:N/2-1)*(Fs/N);

% Find peak frequencies
[peaks, locs] = findpeaks(mag_spectrum, 'MinPeakProminence', max(mag_spectrum)/20);
peak_freqs = freq(locs);

% Sort peaks by magnitude to find dominant frequencies
[sorted_peaks, idx] = sort(peaks, 'descend');
dominant_freqs = peak_freqs(idx);

% Noise analysis
% Estimate noise floor (mean of bottom 10% of magnitudes)
sorted_mags = sort(mag_spectrum);
noise_floor = mean(sorted_mags(1:floor(length(sorted_mags)*0.1)));

% Calculate Signal-to-Noise Ratio (SNR)
signal_power = mean(abs(y).^2);
noise_power = noise_floor^2;
SNR = 10*log10(signal_power/noise_power);

% Plot
figure;
plot(freq, mag_spectrum);
title('Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs / 2]);
grid on;
set(gca, 'YScale', 'log');

% Print audio info
fprintf('Sampling Frequency: %d Hz\n', Fs);
fprintf('Signal Duration: %.2f seconds\n', N/Fs);