function matchBirdSounds(taskFile, referenceFiles)
    chunkSeconds = 0.05; % Chunk length in seconds
    
    [taskY, Fs] = audioread(taskFile);
    
    % Calculate spectrogram parameters
    windowLength = round(chunkSeconds * Fs);
    overlap = round(windowLength / 2);
    nfft = windowLength;
    
    % Get task file spectrogram features
    taskFeatures = extractSpectrogramFeatures(taskY, Fs, windowLength, overlap, nfft);
    
    % Initialize results storage
    similarities = zeros(1, length(referenceFiles));
    
    % Process each reference file
    for i = 1:length(referenceFiles)
        % Load reference file
        [refY, refFs] = audioread(referenceFiles{i});
        
        % Ensure same sampling rate
        if refFs ~= Fs
            refY = resample(refY, Fs, refFs);
        end
        
        % Get reference spectrogram features
        refFeatures = extractSpectrogramFeatures(refY, Fs, windowLength, overlap, nfft);
        
        % Calculate similarity
        similarities(i) = calculateSimilarity2(taskFeatures, refFeatures);
    end
    
    % Find best match
    fprintf('Distance from Bird1: %.2f\n', similarities(1));
    fprintf('Distance from Bird2: %.2f\n', similarities(2));
    % fprintf('Distance from Bird3: %.2f\n', similarities(3));
    % fprintf('Distance from Bird1_1: %.2f\n', similarities(4));
    % fprintf('Distance from Bird2_1: %.2f\n', similarities(5));
    % fprintf('Distance from Bird3_1: %.2f\n', similarities(6));
    % fprintf('Distance from Bird1_2: %.2f\n', similarities(7));
    % fprintf('Distance from Bird2_2: %.2f\n', similarities(8));
    % fprintf('Distance from Bird3_2: %.2f\n', similarities(9));
    
    [maxSim, bestMatch] = min(similarities);
    fprintf('\nBest match is %s with %.2f distance\n', referenceFiles{bestMatch}, maxSim);
    
    %Plot spectrograms side by side
    plotComparison(taskY, audioread(referenceFiles{bestMatch}), Fs, windowLength, overlap, nfft);
end

function features = extractSpectrogramFeatures(y, Fs, windowLength, overlap, nfft)
    % Calculate spectrogram
    [S, ~, ~] = spectrogram(y, hamming(windowLength), overlap, nfft, Fs);
    S = abs(S);  % Magnitude spectrogram
    
    % Normalize spectrogram
    S = S / max(S(:));
    
    % Extract features
    features.S = S;
    features.dominantFreqs = getDominantFrequencies(S, Fs, nfft);
    features.energyProfile = sum(S .* S, 1); % Changed to squared sum
    features.spectralCentroid = calculateSpectralCentroid(S, Fs, nfft);
end

function centroid = calculateSpectralCentroid(S, Fs, nfft)
    freqBins = floor(nfft / 2) + 1;
    frequencies = (0:freqBins-1) * Fs / nfft;
    centroid = sum(frequencies' .* S, 1) ./ sum(S, 1);
end

% Use the findpeaks function instead for this function
function dominantFreqs = getDominantFrequencies(S, Fs, nfft)
    [~, maxFreqIndices] = max(S, [], 1);
    freqBins = floor(nfft / 2) + 1;
    f = (0:freqBins-1) * Fs / nfft;
    dominantFreqs = f(maxFreqIndices);
end

% Removed Spectral Centroid comparison
function similarity = calculateSimilarity2(feat1, feat2)
    % Calculate overall similarity based on multiple features
    
    % 1. Compare dominant frequencies
    freqSim = compareDominantFreqs2(feat1.dominantFreqs, feat2.dominantFreqs);
    
    % 2. Compare energy profiles
    energySim = compareEnergyProfiles2(feat1.energyProfile, feat2.energyProfile);
    
    % 3. Compare spectral centroids
    centroidSim = compareSpectralCentroids(feat1.spectralCentroid, feat2.spectralCentroid);
    
    % Weighted combination of similarities
    weights = [0.00001, 1, 0.00001];  % Adjust weights as needed
    similarity = weights(1)*freqSim + weights(2)*energySim + weights(3)*centroidSim;
end

function sim = compareSpectralCentroids(centroid1, centroid2)
    [dist, ix, iy] = dtw(centroid1, centroid2);

    % figure;
    % subplot(2,1,1);
    % plot(centroid1); hold on; plot(centroid2);
    % title('Original Sequences');
    % legend('Sequence 1', 'Sequence 2');
    % 
    % subplot(2,1,2);
    % plot(1:length(ix), centroid1(ix), 'b', 1:length(iy), centroid2(iy), 'r');
    % title('DTW Aligned Sequences');
    % legend('Aligned Sequence 1', 'Aligned Sequence 2');

    sim = dist;
end

function sim = compareDominantFreqs2(freq1, freq2)
    % Compute DTW distance
    [dist, ix, iy] = dtw(freq1, freq2);

    % figure;
    % subplot(2,1,1);
    % plot(freq1); hold on; plot(freq2);
    % title('Original Sequences');
    % legend('Sequence 1', 'Sequence 2');

    % subplot(2,1,2);
    % plot(1:length(ix), freq1(ix), 'b', 1:length(iy), freq2(iy), 'r');
    % title('DTW Aligned Sequences');
    % legend('Aligned Sequence 1', 'Aligned Sequence 2');

    sim = dist;
end

function sim = compareEnergyProfiles2(energy1, energy2)
    [dist, ix, iy] = dtw(energy1, energy2, "squared");

    % figure;
    % subplot(2,1,1);
    % plot(energy1); hold on; plot(energy2);
    % title('Original Sequences');
    % legend('Sequence 1', 'Sequence 2');

    % subplot(2,1,2);
    % plot(1:length(ix), energy1(ix), 'b', 1:length(iy), energy2(iy), 'r');
    % title('DTW Aligned Sequences');
    % legend('Aligned Sequence 1', 'Aligned Sequence 2');

    sim = dist;
end

function plotComparison(y1, y2, Fs, windowLength, overlap, nfft)
    figure('Name', 'Spectrogram Comparison');
    
    subplot(2, 1, 1);
    spectrogram(y1, hamming(windowLength), overlap, nfft, Fs, 'yaxis');
    colormap('jet');
    title('Task File Spectrogram');
    
    subplot(2, 1, 2);
    spectrogram(y2, hamming(windowLength), overlap, nfft, Fs, 'yaxis');
    colormap('jet');
    title('Best Match Spectrogram');
end

% Define your files
taskFile = 'F4.wav';
% referenceFiles = {'bird1.wav', 'bird2.wav', 'bird3.wav', 'bird1_1.wav', 'bird2_1.wav', 'bird3_1.wav', 'bird1_2.wav', 'bird2_2.wav', 'bird3_2.wav'};
referenceFiles = {'bird1.wav', 'bird3_3.wav'};

% Run the matching
matchBirdSounds(taskFile, referenceFiles);