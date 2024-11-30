function [dominantFreqs, t, f, S] = analyzeAudioFrequencies(audioFile, chunkSeconds, plotResults)
    % ANALYZEAUDIOFREQUENCIES Analyze the dominant frequencies in an audio file
    %   [dominantFreqs, t, f, S] = analyzeAudioFrequencies(audioFile, chunkSeconds, plotResults)
    %
    % Inputs:
    %   audioFile - Path to the audio file (string)
    %   chunkSeconds - Length of each chunk in seconds (default: 0.05)
    %   plotResults - Boolean flag to control plotting (default: false)
    %
    % Outputs:
    %   dominantFreqs - Vector of dominant frequencies over time
    %   t - Time vector
    %   f - Frequency vector
    %   S - Spectrogram matrix
    
    % Set default parameters if not provided
    if nargin < 2
        chunkSeconds = 0.05;
    end
    if nargin < 3
        plotResults = false;
    end
    
    % Read the audio file
    [y, Fs] = audioread(audioFile);
    
    % Convert stereo to mono if necessary
    if size(y, 2) > 1
        y = mean(y, 2);
    end
    
    % Parameters for chunking
    windowLength = round(chunkSeconds * Fs);  % Convert seconds to samples
    overlap = round(windowLength/2);  % Overlap between chunks (50% overlap)
    nfft = windowLength;  % Number of FFT points
    
    % Error checking for window size
    if windowLength > length(y)
        error('Chunk length is too large for the signal');
    end
    
    % Manual implementation parameters
    windowSize = windowLength;
    hopSize = windowLength - overlap;
    
    % Calculate number of frames
    numFrames = floor((length(y) - windowSize)/hopSize) + 1;
    
    % Verify that numFrames is positive
    if numFrames <= 0
        error('Invalid number of frames. Try reducing the chunk length or overlap.');
    end
    
    % Initialize time-frequency matrix
    freqBins = floor(nfft/2) + 1;
    S = zeros(freqBins, numFrames);
    
    % Process each frame
    for i = 1:numFrames
        % Extract frame
        startIdx = (i-1)*hopSize + 1;
        endIdx = min(startIdx + windowSize - 1, length(y));
        frame = y(startIdx:endIdx);
        
        % Zero-pad if necessary
        if length(frame) < windowSize
            frame = [frame; zeros(windowSize - length(frame), 1)];
        end
        
        % Apply window
        frame = frame .* hamming(windowSize);
        
        % Compute FFT
        X = fft(frame, nfft);
        
        % Keep only positive frequencies
        S(:,i) = abs(X(1:freqBins));
    end
    
    % Create time and frequency vectors
    t = (0:numFrames-1)*hopSize/Fs;
    f = (0:freqBins-1)*Fs/nfft;
    
    % Find dominant frequencies
    [~, maxFreqIndices] = max(S, [], 1);
    dominantFreqs = f(maxFreqIndices);
    
    % Apply threshold to ignore silent parts
    % threshold = max(max(S))/10;
    % for i = 1:length(dominantFreqs)
    %     if max(S(:,i)) < threshold
    %         dominantFreqs(i) = NaN;
    %     end
    % end
    
    % Generate plots if requested
    if plotResults
        % Plot spectrogram
        figure;
        subplot(2,1,1);
        imagesc(t, f, 20*log10(S));
        axis xy;
        colormap('jet');
        colorbar;
        title(['Spectrogram (Chunk Length: ' num2str(chunkSeconds*1000) ' ms)']);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        
        % Plot dominant frequencies
        subplot(2,1,2);
        plot(t, dominantFreqs, 'r-', 'LineWidth', 1.5);
        title('Dominant Frequencies Over Time (With Threshold)');
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        grid on;
        ylim([0, Fs/2]);
    end
end

% Basic usage (just get dominant frequencies)
dominantFreqs = analyzeAudioFrequencies('./Reference/bird1.wav');
dominantFreqs2 = analyzeAudioFrequencies('./Task/F2.wav');

dtw(dominantFreqs, dominantFreqs2)