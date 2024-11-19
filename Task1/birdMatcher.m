function matchBirdSounds(taskFile, referenceFiles)
    % Function to compare mono bird sounds and find the best match
    % taskFile: path to the file to be matched
    % referenceFiles: cell array of paths to reference bird sound files

    % Parameters
    chunkSeconds = 0.05;  % Chunk length in seconds
    
    % Load task file
    [taskY, Fs] = audioread(taskFile);
    taskY = ensureMono(taskY);
    
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
        refY = ensureMono(refY);
        
        % Ensure same sampling rate
        if refFs ~= Fs
            refY = resample(refY, Fs, refFs);
        end
        
        % Get reference spectrogram features
        refFeatures = extractSpectrogramFeatures(refY, Fs, windowLength, overlap, nfft);
        
        % Calculate similarity
        similarities(i) = calculateSimilarity(taskFeatures, refFeatures);
        
        % Display individual similarity
        fprintf('Similarity with %s: %.2f%%\n', referenceFiles{i}, similarities(i) * 100);
    end
    
    % Find best match
    [maxSim, bestMatch] = max(similarities);
    fprintf('\nBest match is %s with %.2f%% similarity\n', ...
        referenceFiles{bestMatch}, maxSim * 100);
    
    % Optional: Plot spectrograms side by side
    plotComparison(taskY, audioread(referenceFiles{bestMatch}), Fs, windowLength, overlap, nfft);
end

function y = ensureMono(y)
    % Ensure the signal is mono by averaging channels if stereo
    if size(y, 2) > 1
        y = mean(y, 2);
    end
end

function features = extractSpectrogramFeatures(y, Fs, windowLength, overlap, nfft)
    % Calculate spectrogram
    hopSize = windowLength - overlap;
    [S, ~, ~] = spectrogram(y, hamming(windowLength), overlap, nfft, Fs);
    S = abs(S);  % Magnitude spectrogram
    
    % Normalize spectrogram
    S = S / max(S(:));
    
    % Extract features
    features.S = S;
    features.dominantFreqs = getDominantFrequencies(S, Fs, nfft);
    features.energyProfile = sum(S, 1);
    features.spectralCentroid = calculateSpectralCentroid(S, Fs, nfft);
end

function dominantFreqs = getDominantFrequencies(S, Fs, nfft)
    [~, maxFreqIndices] = max(S, [], 1);
    freqBins = floor(nfft / 2) + 1;
    f = (0:freqBins-1) * Fs / nfft;
    dominantFreqs = f(maxFreqIndices);
    
    % Apply threshold
    threshold = max(S(:)) / 10;
    dominantFreqs(max(S, [], 1) < threshold) = NaN;
end

function centroid = calculateSpectralCentroid(S, Fs, nfft)
    freqBins = floor(nfft / 2) + 1;
    frequencies = (0:freqBins-1) * Fs / nfft;
    centroid = sum(frequencies' .* S, 1) ./ sum(S, 1);
end

function similarity = calculateSimilarity(feat1, feat2)
    % Calculate overall similarity based on multiple features
    
    % 1. Compare dominant frequencies
    freqSim = compareDominantFreqs(feat1.dominantFreqs, feat2.dominantFreqs);
    
    % 2. Compare energy profiles
    energySim = compareEnergyProfiles(feat1.energyProfile, feat2.energyProfile);
    
    % 3. Compare spectral centroids
    centroidSim = compareSpectralCentroids(feat1.spectralCentroid, feat2.spectralCentroid);
    
    % Weighted combination of similarities
    weights = [0.4, 0.3, 0.3];  % Adjust weights as needed
    similarity = weights(1)*freqSim + weights(2)*energySim + weights(3)*centroidSim;
end

function sim = compareDominantFreqs(freq1, freq2)
    freq1 = freq1(~isnan(freq1));
    freq2 = freq2(~isnan(freq2));
    
    len = min(length(freq1), length(freq2));
    if len > 1
        freq1 = interp1(1:length(freq1), freq1, linspace(1, length(freq1), len), 'linear', 'extrap');
        freq2 = interp1(1:length(freq2), freq2, linspace(1, length(freq2), len), 'linear', 'extrap');
    end
    
    diff = abs(freq1 - freq2);
    sim = 1 - mean(diff) / max([freq1(:); freq2(:)]);
end

function sim = compareEnergyProfiles(energy1, energy2)
    len = min(length(energy1), length(energy2));
    energy1 = interp1(1:length(energy1), energy1, linspace(1, length(energy1), len), 'linear', 'extrap');
    energy2 = interp1(1:length(energy2), energy2, linspace(1, length(energy2), len), 'linear', 'extrap');
    
    energy1 = energy1 / max(energy1);
    energy2 = energy2 / max(energy2);
    
    sim = max(xcorr(energy1, energy2, 'normalized'));
end

function sim = compareSpectralCentroids(centroid1, centroid2)
    len = min(length(centroid1), length(centroid2));
    centroid1 = interp1(1:length(centroid1), centroid1, linspace(1, length(centroid1), len), 'linear', 'extrap');
    centroid2 = interp1(1:length(centroid2), centroid2, linspace(1, length(centroid2), len), 'linear', 'extrap');
    
    diff = abs(centroid1 - centroid2);
    sim = 1 - mean(diff) / max([centroid1(:); centroid2(:)]);
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
taskFile = 'Asine.wav';
referenceFiles = {'bird1.wav', 'bird2.wav', 'bird3.wav'};

% Run the matching
matchBirdSounds(taskFile, referenceFiles);