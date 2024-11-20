% Read the audio file
% [audio, fs] = audioread(fullfile('..','References and Tasks', 'Project_LouderWordsDetection', 'audios', '1.wav'));  % Replace with your wav filename

% % Read timestamp data from text file
% opts = detectImportOptions('1.txt', 'Delimiter', '\t');
% opts.VariableNames = {'Word', 'StartTime', 'EndTime', 'Expected'};
% opts = setvartype(opts, {'Word', 'StartTime', 'EndTime', 'Expected'}, {'char', 'double', 'double', 'double'});
% data = readtable('1.txt', opts);

% words = data.Word;
% startTimes = data.StartTime;
% endTimes = data.EndTime;
% expected = data.Expected;
[audio, fs] = audioread('8.mp3'); 

[words, startTimes, endTimes, isLouder] = readTimeFile('8.txt');

% Take the FFT of the audio signal
n = length(audio);  % Number of samples
f = (0:n-1)*(fs/n);  % Frequency range
y = fft(audio);  % Compute the FFT

% Plot the magnitude of the FFT
figure;
plot(f, abs(y));
title('Magnitude of FFT of Audio Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 fs/2]);  % Plot up to the Nyquist frequency
grid on;

% Print to verify the data
for i = 1:length(words)
    fprintf('%s\t\t%f\t%f\t%d\n', words{i}, startTimes(i), endTimes(i), isLouder(i));
end
fprintf('\n');

% Process each word
for i = 1:length(words)
    % Convert time to samples
    startSample = round(startTimes(i) * fs);
    endSample = round(endTimes(i) * fs);
    
    % Extract word segment    
    wordSegment = audio(startSample:endSample);
    duration = endTimes(i) - startTimes(i);
    
    % Calculate RMS energy
    energy = sqrt(mean(wordSegment.^2));
    peakAmplitude = max(abs(wordSegment));
    normalizedEnergy = energy / duration;
    % Perform STFT
    [S, F, T] = stft(wordSegment, fs, 'Window', hamming(256), 'OverlapLength', 128, 'FFTLength', 512);

    % Calculate energy in the frequency band of interest (40-500Hz)
    freqBand = (F >= 40 & F <=500);
    bandEnergy = sum(abs(S(freqBand, :)).^2, 'all');

    normalizedBandEnergy = bandEnergy / duration;
    % NOTE: 'bandEnergy' is the energy in the frequency domain but the 'energy' is in time domain.

    % Determine if loud (you may need to adjust threshold)
    threshold = 0.12;  % Adjust based on your audio
    isLoud = energy > threshold;
    
    % Print result
    fprintf('Word: %s \t Peak Amplitude: %.4f \t Energy: %.4f \t normalisedEnergy: %.4f \t Band_Energy: %.4f \t normalisedBand_Energy: %.4f \t Is Loud: %d\n', words{i}, peakAmplitude, energy, normalizedEnergy, bandEnergy, normalizedBandEnergy, isLoud);
end

% The above code gives the audio characteristics like peak amplitude, energy, normalized energy, band energy, normalized band energy to determine if the word is loud or not.
% Assumption: The threshold value is set to 0.12 currently (can adjust it later)
% the band for human voice is taken to be from 40 to 500 Hz (got this idea by plotting by FFT of the audio signal)
% Observation: if amplitude > 0.6, sure shot 100% loud word (based on the data given)
% Check the word told in the 6th txt, parameters shows that it shows be loud