function energyOverTime = energyPlot(audio, fs)
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

function [wordTimestamps] = detectWordTimestamps(audioFilePath)
    % Word Timestamp Detection in Audio Signal
    % Input: Path to audio file
    % Output: Array of word start and end timestamps
    
    % Read the audio file
    [audioSignal, sampleRate] = audioread(audioFilePath);
    
    % Convert to mono if stereo
    if size(audioSignal, 2) > 1
        audioSignal = mean(audioSignal, 2);
    end

    energy = energyPlot(audioSignal, sampleRate);

    % Threshold for voice activity (adjust as needed)
    energyThreshold = 0.5;
    
    % Detect speech regions
    isSpeech = energy > energyThreshold;
    
    % Find speech segment boundaries
    speechOnsets = find(diff([0, isSpeech]) == 1);
    speechOffsets = find(diff([isSpeech, 0]) == -1);
    
    % Convert sample indices to timestamps
    wordTimestamps = [];
    for i = 1:length(speechOnsets)
        onsetTime = (speechOnsets(i)) / sampleRate;
        offsetTime = (speechOffsets(i)) / sampleRate;
        
        % Optional: Apply additional filtering for short/long segments
        segmentDuration = offsetTime - onsetTime;
        if segmentDuration > 0.01 && segmentDuration < 4.0
            wordTimestamps(end+1, 1:2) = [onsetTime, offsetTime];
        end
    end
end

% Example usage
audioFilePath = '2.wav';
timestamps = detectWordTimestamps(audioFilePath);
disp('Word Timestamps:');
disp(timestamps);