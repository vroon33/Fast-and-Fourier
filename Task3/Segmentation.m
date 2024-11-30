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
    
    % Normalize the audio signal
    audioSignal = audioSignal / max(abs(audioSignal));
    
    % Energy-based Voice Activity Detection (VAD)
    windowSize = round(0.025 * sampleRate);  % 25 ms window
    hopLength = round(0.01 * sampleRate);    % 10 ms hop
    
    % Calculate short-time energy
    energy = zeros(1, floor((length(audioSignal) - windowSize) / hopLength) + 1);
    for i = 1:length(energy)
        startIdx = 1 + (i-1) * hopLength;
        endIdx = startIdx + windowSize - 1;
        energy(i) = sum(audioSignal(startIdx:endIdx).^2);
    end
    
    % Normalize energy
    energy = (energy - min(energy)) / (max(energy) - min(energy));
    
    % Threshold for voice activity (adjust as needed)
    energyThreshold = 0.001;
    
    % Detect speech regions
    isSpeech = energy > energyThreshold;
    
    % Find speech segment boundaries
    speechOnsets = find(diff([0, isSpeech]) == 1);
    speechOffsets = find(diff([isSpeech, 0]) == -1);
    
    % Convert sample indices to timestamps
    wordTimestamps = [];
    for i = 1:length(speechOnsets)
        onsetTime = (speechOnsets(i) * hopLength) / sampleRate;
        offsetTime = (speechOffsets(i) * hopLength) / sampleRate;
        
        % Optional: Apply additional filtering for short/long segments
        segmentDuration = offsetTime - onsetTime;
        if segmentDuration > 0.1 && segmentDuration < 2.0
            wordTimestamps(end+1, 1:2) = [onsetTime, offsetTime];
        end
    end
     
    % Optional: Additional refinement using spectral subtraction
    % This is a basic implementation and might need tuning
    function refinedTimestamps = refineTimestamps(timestamps, audioSignal, sampleRate)
        refinedTimestamps = timestamps;
        
        % You could add more sophisticated word boundary detection here
        % Examples:
        % - Zero-crossing rate analysis
        % - Spectral flux
        % - More advanced silence detection
    end
    
    % Optional post-processing
    wordTimestamps = refineTimestamps(wordTimestamps, audioSignal, sampleRate);
end

% Example usage
audioFilePath = '1.wav';
timestamps = detectWordTimestamps(audioFilePath);
disp('Word Timestamps:');
disp(timestamps);