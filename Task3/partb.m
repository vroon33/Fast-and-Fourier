clear all; close all;
function [energy,samplerate] = energycal(y, fs)
    % Parameters for analysis
    samplerate = 0.004;
    windowSize = round(samplerate * fs);  % 10 ms window
    
    % Overlap between windows
    overlap = 0;
    
    % Initialize energy array
    energy = zeros(1, ceil(length(y) / (windowSize - overlap)));
    
    % Calculate energy for each window
    for i = 1:length(energy)
        % Calculate start and end of current window
        startIdx = 1 + (i-1) * (windowSize - overlap);
        endIdx = min(startIdx + windowSize - 1, length(y));
        
        % Extract window
        window = y(startIdx:endIdx);
        
        % Calculate instantaneous energy (squared amplitude)
        energy(i) = sum(window.^2);
    end
    
    % Normalize energy
    energy = energy / max(energy);
end

% Read the audio file
    [y, fs] = audioread('2.wav');
    % y = noisefilter(y, fs);
    
    % Convert stereo to mono if necessary
    if size(y, 2) > 1
        y = mean(y, 2);
    end
    
    % Noise filtering
    % y = noisefilter(y, fs); 

    [energy,samplerate] = energycal(y, fs);

    threshold = 0.1;
    % Find segments with high energy (emphasized parts)
    isSpeech = energy > threshold * mean(energy);

    % Find speech segment boundaries
    speechOnsets = find(diff([0, isSpeech]) == 1);
    speechOffsets = find(diff([isSpeech, 0]) == -1);
    
    % Convert sample indices to timestamps
    wordTimestamps = [];
    onsetTime = speechOnsets .* samplerate;
    offsetTime = speechOffsets .* samplerate;

    % Old Approach
    % i = 1;
    % while i <= length(speechOnsets)
    %     %Applying additional filtering for short/long segments
    %     segmentDuration = offsetTime(i) - onsetTime(i);
        
    %     if segmentDuration > 0.08 && segmentDuration < 4.0
    %         wordTimestamps(end+1, 1:2) = [onsetTime(i), offsetTime(i)];
    %         i= i+1;
    %     else 
    %         i = i+1;
    %     end
    % end

    % middle approach
    i = 1;
    while i <= length(speechOnsets)
        %Applying additional filtering for short/long segments
        segmentDuration = offsetTime(i) - onsetTime(i);
        
        if segmentDuration > 0.045 && segmentDuration < 4.0
            if i < length(speechOnsets) && (onsetTime(i+1) - offsetTime(i) < 0.02)
            % Combine segments
            wordTimestamps(end+1, 1:2) = [onsetTime(i), offsetTime(i+1)];
            i=i+2;
            else
            wordTimestamps(end+1, 1:2) = [onsetTime(i), offsetTime(i)];
            i= i+1;
            end
        else 
            i = i+1;
        end
    end

    
    % New Approach
    % i = 1;
    % while i <= length(speechOnsets)
    %     segmentDuration = offsetTime(i) - onsetTime(i);
        
    %     % Check if there's a next segment and if it's close to current segment
    %     if i < length(speechOnsets) && (onsetTime(i+1) - offsetTime(i)) < 0.08
    %         % Combine segments
    %         wordTimestamps(end+1, 1:2) = [onsetTime(i), offsetTime(i+1)];
    %         i = i + 2; % Skip next segment since we've combined it
    %     elseif segmentDuration > 0.08 && segmentDuration < 4.0
    %         % Regular segment processing
    %         wordTimestamps(end+1, 1:2) = [onsetTime(i), offsetTime(i)];
    %         i = i + 1;
    %     else
    %         i = i + 1;
    %     end
    % end

    emphasized = zeros(size(energy));
    for i = 1:size(wordTimestamps, 1)
        startIdx = round(wordTimestamps(i,1)/samplerate);
        endIdx = round(wordTimestamps(i,2)/samplerate);
        emphasized(startIdx:endIdx) = 1;
    end

    % Create time axis
    timeAxis = ((0:length(energy)-1) * (round(samplerate*fs))) / fs;
    % Plot results
    figure;
    
    % Plot waveform
    % subplot(3,1,1);
    % plot((1:length(y))/fs, y);
    % title('Audio Waveform');
    % xlabel('Time (s)');
    % ylabel('Amplitude');
    
    % Plot energy
    subplot(3,1,1);
    plot(timeAxis, energy);
    title('Energy vs Time');
    xlabel('Time (seconds)');
    ylabel('Energy (Squared Amplitude)');
    grid on;
    
    % plot isSpeech
    subplot(3,1,2)
    plot(timeAxis, isSpeech);
    title('Speech Detection');
    xlabel('Time (s)');
    ylabel('Speech Detection');
    ylim([-0.1 1.1]);

    % Plot emphasized segments
    subplot(3,1,3);
    plot(timeAxis, emphasized);
    title('Emphasized Segments (1 = emphasized)');
    xlabel('Time (s)');
    ylabel('Emphasis Detection');
    ylim([-0.1 1.1]);
    
    % % Print time stamps of emphasized segments
    % emphasisChanges = diff([0; emphasized; 0]);
    % startPoints = find(emphasisChanges == 1);
    % endPoints = find(emphasisChanges == -1) - 1;
    
    % fprintf('\nEmphasized segments found at:\n');
    % for i = 1:length(startPoints)
    %     startTime = timeAxis(startPoints(i));
    %     endTime = timeAxis(endPoints(i));
    %     fprintf('%.2f seconds to %.2f seconds\n', startTime, endTime);
    % end
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% CODE FROM PART A

% Arrays to store metrics for each word
peakAmplitudes = zeros(1, length(wordTimestamps));
energies = zeros(1, length(wordTimestamps));
normalizedEnergies = zeros(1, length(wordTimestamps));
bandEnergies = zeros(1, length(wordTimestamps));
normalizedBandEnergies = zeros(1, length(wordTimestamps));

% Process each word - first pass to collect metrics
for i = 1:length(wordTimestamps)
    startSample = round(wordTimestamps(i,1) * fs);
    endSample = round(wordTimestamps(i,2) * fs);
    wordSegment = y(startSample:endSample);
    duration = wordTimestamps(i,2) - wordTimestamps(i,1);
    
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
for i = 1:length(wordTimestamps)
    % Condition 3: Check if peak amplitude is greater than mean
    isPeakHigher = peakAmplitudes(i) > meanPeakAmplitude;
    
    % Condition 1: Check if energy is greater than 2 * mean energy
    isEnergyHigher = energies(i) > 2 * meanEnergy;

    % Condition 4: Check if normalised band energy is greater than 2*mean normalised band energy
    isBandEnergyHigher = normalizedBandEnergies(i) > 2 * meanNormalizedBandEnergy;
    
    % Apply the conditions in order
    isLoud = 0;
    if isPeakHigher  % First check condition 3
        if isBandEnergyHigher  % Then check condition 4
            isLoud = 1;
        elseif isEnergyHigher  % Then check condition 1
            isLoud = 1;
        end
    end
    fprintf('From t=%.2f to %.2f \t Peak Amplitude: %.4f \t Energy: %.4f \t Is Loud: %d\n', wordTimestamps(i,1), wordTimestamps(i,2), peakAmplitudes(i), energies(i), isLoud);
end