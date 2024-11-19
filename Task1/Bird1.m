% %% Script to find the FFT of bird1.wav

% % Read the audio file
% [y, Fs] = audioread('bird1.wav');

% % If stereo, convert to mono by averaging channels
% if size(y, 2) > 1
%     y = mean(y, 2);
% end

% % Get the length of the signal
% N = length(y);

% % Compute the FFT
% Y = fft(y);

% % Compute the two-sided spectrum
% P2 = abs(Y/N);

% % Compute the single-sided spectrum
% P1 = P2(1:floor(N/2)+1);
% P1(2:end-1) = 2*P1(2:end-1);

% % Define the frequency domain
% f = Fs*(0:(N/2))/N;

% % Create the frequency spectrum plot
% figure;
% plot(f, P1);
% title('Single-Sided Amplitude Spectrum');
% xlabel('Frequency (Hz)');
% ylabel('|P1(f)|');
% grid on;

% % Adjust the x-axis to show frequencies up to Fs/2
% xlim([0 Fs/2]);

% % Add logarithmic scale for better visualization of amplitude
% set(gca, 'YScale', 'log');

% % Optional: Add information about the signal
% fprintf('Sampling Frequency: %d Hz\n', Fs);
% fprintf('Signal Duration: %.2f seconds\n', N/Fs);

% Read the audio file
[y, Fs] = audioread('Bird1.wav');

% Check if the audio is mono or stereo
[~, numChannels] = size(y);
isStereo = (numChannels > 1);

% Create the figure with appropriate number of subplots
if isStereo
    figure('Position', [100 100 1500 500]);
    numPlots = 3;  % Three plots for stereo (mono + left + right)
else
    figure('Position', [100 100 500 500]);
    numPlots = 1;  % One plot for mono
end

% Get signal length
N = length(y);
f = Fs*(0:(N/2))/N;  % Frequency vector

if isStereo
    % First plot: Combined/Mono signal
    subplot(1,3,1);
    y_mono = mean(y, 2);  % Convert to mono by averaging channels
else
    y_mono = y;  % Already mono
end

% Compute FFT for mono signal
Y_mono = fft(y_mono);
P2_mono = abs(Y_mono/N);
P1_mono = P2_mono(1:floor(N/2)+1);
P1_mono(2:end-1) = 2*P1_mono(2:end-1);

% Plot mono spectrum
if isStereo
    subplot(1,3,1);
else
    subplot(1,1,1);
end
plot(f, P1_mono);
title('Combined (Mono) Spectrum');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
grid on;
xlim([0 Fs/2]);
set(gca, 'YScale', 'log');

if isStereo
    % Process and plot left channel (channel 1)
    subplot(1,3,2);
    Y_left = fft(y(:,1));
    P2_left = abs(Y_left/N);
    P1_left = P2_left(1:floor(N/2)+1);
    P1_left(2:end-1) = 2*P1_left(2:end-1);

    plot(f, P1_left);
    title('Left Channel Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');
    grid on;
    xlim([0 Fs/2]);
    set(gca, 'YScale', 'log');

    % Process and plot right channel (channel 2)
    subplot(1,3,3);
    Y_right = fft(y(:,2));
    P2_right = abs(Y_right/N);
    P1_right = P2_right(1:floor(N/2)+1);
    P1_right(2:end-1) = 2*P1_right(2:end-1);

    plot(f, P1_right);
    title('Right Channel Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');
    grid on;
    xlim([0 Fs/2]);
    set(gca, 'YScale', 'log');
end

% Add appropriate main title
if isStereo
    sgtitle('Frequency Spectrum Analysis - Complete View (Stereo)');
else
    sgtitle('Frequency Spectrum Analysis - Mono Signal');
end

% Print audio information
fprintf('\nAudio File Information:\n');
fprintf('Sampling Frequency: %d Hz\n', Fs);
fprintf('Signal Duration: %.2f seconds\n', N/Fs);
fprintf('Number of Channels: %d\n', numChannels);

% Optional: Add peak detection
[peaks_mono, locs_mono] = findpeaks(P1_mono, f, 'MinPeakProminence', max(P1_mono)/10);
fprintf('\nDominant Frequencies (Combined Signal):\n');
for i = 1:min(5, length(locs_mono))
    fprintf('%.1f Hz (Magnitude: %.2e)\n', locs_mono(i), peaks_mono(i));
end