% Final Super Ultimate Promax Spectacular Code for Task 3a
% Task 3a Completed!
% ab to sahi mai complete ho gaya!
clear all; close all;

function stereoToMono(audio, Fs)
    % Convert stereo to mono if necessary
    if size(audio, 2) > 1
        audio = mean(audio, 2);
    end
end

function plotfft(audio,fs)
    % Take the FFT of the audio signal
    n = length(audio);  % Number of samples
    f = (0:n-1)*(fs/n);  % Frequency range
    y = fft(audio);  % Compute the FFT

    % Plot the magnitude of the FFT
    % figure;
    plot(f, abs(y));
    title('Magnitude of FFT of Audio Signal');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 fs/2]);  % Plot up to the Nyquist frequency
    grid on;
end

function plotsig(audio, fs)
    % Plot the audio signal
    t = (0:length(audio)-1) / fs;  % Time vector
    plot(t, audio);
    % title('Audio Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

function filteredAudio = filter(audio,fs)
    % Plays the audio signal
    sound(audio, fs);

    % Perform FFT on the audio signal
    Y = fft(audio);

    % Get the length of the audio and calculate frequency vector
    N = length(audio);
    freq = (0:N-1) * (fs / N);

    % Create a mask for the desired frequency range
    mask = (freq <= 7000) & (freq >= 100);
    % For negative frequencies (second half of the FFT)
    mask = mask | (freq >= fs - 7000 & freq <= fs - 100);

    figure;
    % Plot the frequency domain of the original audio
    subplot(2,1,1);
    plot(freq, abs(Y));
    title('Frequency Domain of Original Audio');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 fs/2]);  % Plot up to the Nyquist frequency
    grid on;

    % Zero out frequencies outside the desired range
    Y(~mask) = 0;

    % Plot the frequency domain of Y
    subplot(2,1,2);
    plot(freq, abs(fft(filteredAudio)));
    title('Frequency Domain of Filtered Audio');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 fs/2]);  % Plot up to the Nyquist frequency
    grid on;

    % Inverse FFT to get back to time domain
    filteredAudio = real(ifft(Y));

    % Normalize the audio to prevent clipping
    filteredAudio = filteredAudio / max(abs(filteredAudio));

    % Play the filtered audio
    sound(filteredAudio, fs);
end

function energyPlot(audio, fs)
    % Window size (adjust as needed)
    windowSize = round(0.01 * fs);  % 10 ms window
    
    % Overlap between windows
    % overlap = round(windowSize / 2);
    overlap = 0;
    
    % Initialize energy array
    energyOverTime = zeros(1, ceil(length(audio) / (windowSize - overlap)));
    
    % Calculate energy for each window
    for i = 1:length(energyOverTime)
        % Calculate start and end of current window
        startIdx = 1 + (i-1) * (windowSize - overlap);
        endIdx = min(startIdx + windowSize - 1, length(audio));
        
        % Extract window
        window = audio(startIdx:endIdx);
        
        % Calculate instantaneous energy (squared amplitude)
        energyOverTime(i) = sum(window.^2);
    end
    
    % Create time axis
    timeAxis = ((0:length(energyOverTime)-1) * (windowSize - overlap)) / fs;
    
    % Plot energy over time
    figure;
    plot(timeAxis, energyOverTime);
    title('Energy vs Time');
    xlabel('Time (seconds)');
    ylabel('Energy (Squared Amplitude)');
    grid on;
end

[audio, fs] = audioread('9.mp3'); 
[words, startTimes, endTimes, isLouder] = readTimeFile('9.txt');

stereoToMono(audio, fs);

% energyPlot(audio, fs);
% % Plot the original audio signal
% figure;
% subplot(2, 1, 1);
% plotsig(audio, fs);
% title('Original Audio Signal');

% % Plot the filtered audio signal
% subplot(2, 1, 2);
% plotsig(filteredAudio, fs);
% title('Filtered Audio Signal');

% figure;
% subplot(2, 1, 1);
% plotfft(audio, fs);
% subplot(2, 1, 2);
% plotfft(filteredAudio, fs);

% Arrays to store metrics for each word
peakAmplitudes = zeros(1, length(words));
energies = zeros(1, length(words));
normalizedEnergies = zeros(1, length(words));
bandEnergies = zeros(1, length(words));
normalizedBandEnergies = zeros(1, length(words));

% Process each word - first pass to collect metrics
for i = 1:length(words)
    startSample = round(startTimes(i) * fs);
    endSample = round(endTimes(i) * fs);
    wordSegment = audio(startSample:endSample);
    duration = endTimes(i) - startTimes(i);
    
    % Calculate metrics
    energies(i) = sum(wordSegment.^2);
    peakAmplitudes(i) = max(abs(wordSegment));
    normalizedEnergies(i) = energies(i) / duration;

    [S, F, T] = stft(wordSegment, fs, 'Window', hamming(256), 'OverlapLength', 128, 'FFTLength', 512);
    freqBand = (F >= 100 & F <= 7000);
    bandEnergies(i) = sum(abs(S(freqBand, :)).^2, 'all');
    normalizedBandEnergies(i) = bandEnergies(i) / duration;
end

% Calculate means
meanPeakAmplitude = mean(peakAmplitudes);
meanEnergy = mean(energies);
meanNormalizedEnergy = mean(normalizedEnergies);
meanNormalizedBandEnergy = mean(normalizedBandEnergies);

fprintf('\nMean Values:\n');
fprintf('Mean Peak Amplitude: %.5f\n', meanPeakAmplitude);
fprintf('Mean Energy: %.5f\n', meanEnergy);
fprintf('Mean normalizedEnergy: %.5f\n', meanNormalizedEnergy);
fprintf('Mean normalizedBand_Energy: %.5f\n', meanNormalizedBandEnergy);

% Second pass - determine loudness using new conditions
fprintf('\nWord Analysis:\n');
for i = 1:length(words)
    % Condition 3: Check if peak amplitude is greater than mean
    isPeakHigher = peakAmplitudes(i) > meanPeakAmplitude;
    
    % Condition 1: Check if energy is greater than 2 * mean energy
    isEnergyHigher = energies(i) > 2 * meanEnergy;

    % Condition 4: Check if normalised band energy is greater than 2*mean normalised band energy
    isBandEnergyHigher = normalizedBandEnergies(i) > 2 * meanNormalizedBandEnergy;
    
    % Apply the conditions in order
    isLoud = 0;
    if isPeakHigher  % First check condition 3
        if isBandEnergyHigher || isEnergyHigher  % Then check condition 4 and 1
            isLoud = 1;
        end
    end
    fprintf('Word: %s \t Peak Amplitude: %.4f \t Energy: %.4f \t Is Loud: %d\n', words{i}, peakAmplitudes(i), energies(i), isLoud);
end

