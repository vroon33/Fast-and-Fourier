% Read the audio file
[y, Fs] = audioread('F5.wav');

% Convert stereo to mono if necessary
if size(y, 2) > 1
    y = mean(y, 2);
end

% Parameters for chunking
chunkSeconds = 0.05;  % Specify chunk length in seconds (modify this value)
windowLength = round(chunkSeconds * Fs);  % Convert seconds to samples
overlap = round(windowLength/2);  % Overlap between chunks (50% overlap)
nfft = windowLength;  % Number of FFT points

% Error checking for window size
if windowLength > length(y)
    error('Chunk length is too large for the signal');
end

% Create spectrogram using built-in function
figure;
spectrogram(y, hamming(windowLength), overlap, nfft, Fs, 'yaxis');
colormap('jet');
colorbar;

% Customize the plot
title(['Spectrogram of Bird Sound (Chunk Length: ' num2str(chunkSeconds*1000) ' ms)']);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Manual implementation
windowSize = windowLength;
hopSize = windowLength - overlap;

% Calculate number of frames
numFrames = floor((length(y) - windowSize)/hopSize) + 1;

% Verify that numFrames is positive
if numFrames <= 0
    error('Invalid number of frames. Try reducing the chunk length or overlap.');
end

% Display frame information
disp(['Number of frames: ' num2str(numFrames)]);

% Initialize time-frequency matrix
freqBins = floor(nfft/2) + 1;  % Ensure integer number of frequency bins
S = zeros(freqBins, numFrames);

% Process each frame
for i = 1:numFrames
    % Extract frame
    startIdx = (i-1)*hopSize + 1;
    endIdx = min(startIdx + windowSize - 1, length(y));  % Ensure we don't exceed signal length
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

% Plot manual spectrogram
figure;
imagesc(t, f, 20*log10(S));
axis xy;
colormap('jet');
colorbar;
title(['Manual Spectrogram of Bird Sound (Chunk Length: ' num2str(chunkSeconds*1000) ' ms)']);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Display some useful information
disp(['Sampling Frequency: ' num2str(Fs) ' Hz']);
disp(['Chunk Length: ' num2str(chunkSeconds*1000) ' ms']);
disp(['Samples per chunk: ' num2str(windowLength)]);
disp(['Frequency Resolution: ' num2str(Fs/windowLength) ' Hz']);
disp(['Time Resolution: ' num2str(hopSize/Fs*1000) ' ms']);

% Optional: Play the sound
% sound(y, Fs);