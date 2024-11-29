% Load Reference files (bird1.wav, bird2.wav, bird3.wav)
[bird1, fs] = audioread('bird1.wav');
[bird2, ~] = audioread('bird2.wav');
[bird3, ~] = audioread('bird3.wav');

% Load Task file (input dynamically)
taskFile = 'F7.wav'; % Replace with desired task file name
[taskAudio, ~] = audioread(taskFile);

% Segment duration in seconds
segmentDuration = 0.01; % 100 ms
segmentLength = round(segmentDuration * fs);

% Analyze all files
bird1_features = extractDominantFrequencies(bird1, fs, segmentLength);
bird2_features = extractDominantFrequencies(bird2, fs, segmentLength);
bird3_features = extractDominantFrequencies(bird3, fs, segmentLength);

% Extract features for the Task file
task_features = extractDominantFrequencies(taskAudio, fs, segmentLength);

% Normalize and compare features
bird1_prob = calculateSimilarity(bird1_features, task_features);
bird2_prob = calculateSimilarity(bird2_features, task_features);
bird3_prob = calculateSimilarity(bird3_features, task_features);

% Find the best match using total probabilities
[~, bestMatch] = max([bird1_prob, bird2_prob, bird3_prob]);

% Output the result
fprintf('The Task file (%s) matches best with Bird Species %d\n', taskFile, bestMatch);

% Function to extract dominant frequencies
function features = extractDominantFrequencies(audio, fs, segmentLength)
    % Preallocate storage for features
    numSegments = floor(length(audio) / segmentLength);
    features = zeros(numSegments, 4); % Nx4 matrix for 4 peak frequencies

    for i = 1:numSegments
        % Extract segment
        segment = audio((i-1)*segmentLength + 1 : i*segmentLength);

        % Apply FFT and filter
        fftData = abs(fft(segment));
        fftData = fftData(1:floor(end/2)); % Single-sided spectrum
        freqs = (0:length(fftData)-1) * (fs / length(segment));
        
        % Remove very low and very high frequencies (5 Hz to 15 kHz)
        validIdx = freqs > 5 & freqs < 15000;
        fftData = fftData(validIdx);
        freqs = freqs(validIdx);

        % Normalize
        fftData = fftData / max(fftData);

        % Find top 4 peak frequencies
        [~, peakIdx] = sort(fftData, 'descend');
        peakFreqs = freqs(peakIdx(1:4));

        % Store in the matrix
        features(i, :) = peakFreqs';
    end
end

% Function to calculate similarity (total probability)
function prob = calculateSimilarity(referenceFeatures, taskFeatures)
    % Normalize features
    refNorm = normalize(referenceFeatures, 'range', [0, 1]);
    taskNorm = normalize(taskFeatures, 'range', [0, 1]);

    % Compare using product of probabilities
    prob = 1;
    for i = 1:size(taskNorm, 2)
        % Calculate element-wise similarity
        distances = abs(refNorm(:, i) - taskNorm(:, i));
        minDist = min(distances);
        prob = prob * (1 - minDist); % Convert distance to probability
    end
end
