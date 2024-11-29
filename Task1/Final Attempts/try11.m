% List the reference audio files
ref_files = {'bird1.wav', 'bird2.wav', 'bird3.wav'};

% Load the task audio file
task_file = 'F1.wav';
[task_audio, task_fs] = audioread(task_file);

% Threshold for the score percentile (fraction of max possible value)
threshold = 0.9;

% Initialize variables to track the best match
best_match = '';
best_score = -inf;
method_used = '';

% Initialize score arrays
dominant_freq_scores = zeros(1, length(ref_files));
spcc_scores = zeros(1, length(ref_files));
tdcc_scores = zeros(1, length(ref_files));

% Stage 1: Dominant Frequency Matching + Spectrogram Cross-Correlation (SPCC)
% fprintf('Stage 1: Dominant Frequencies + Spectrogram Cross-Correlation (SPCC)\n'); % Commented out
for j = 1:length(ref_files)
    ref_file = ref_files{j};
    [ref_audio, ref_fs] = audioread(ref_file);
    
    % Compute spectrograms
    [task_spectrogram, task_frequencies, ~] = spectrogram(task_audio, hamming(256), 128, 256, task_fs);
    [ref_spectrogram, ref_frequencies, ~] = spectrogram(ref_audio, hamming(256), 128, 256, ref_fs);
    
    % Compute dominant frequencies and their order
    [task_dominant_freqs, task_dominant_order] = find_dominant_frequencies(task_spectrogram, task_frequencies, 4);
    [ref_dominant_freqs, ref_dominant_order] = find_dominant_frequencies(ref_spectrogram, ref_frequencies, 4);
    
    % Compute similarity score based on dominant frequencies
    dominant_freq_scores(j) = compute_similarity_score(ref_dominant_freqs, ref_dominant_order, task_dominant_freqs, task_dominant_order);
    
    % Compute spectrogram cross-correlation
    [spcc, ~] = xcorr(task_spectrogram(:), ref_spectrogram(:));
    spcc_scores(j) = max(spcc);
    
    % fprintf('    Match with %s: Dominant Freq Score: %.4f, SPCC Score: %.4f\n', ref_file, dominant_freq_scores(j), spcc_scores(j)); % Commented out
end

% Stage 2: Time-Domain Cross-Correlation (TDCC)
% fprintf('\nStage 2: Time-Domain Cross-Correlation (TDCC)\n'); % Commented out
for j = 1:length(ref_files)
    ref_file = ref_files{j};
    [ref_audio, ~] = audioread(ref_file);
    
    % Compute TDCC
    [tdcc, ~] = xcorr(task_audio, ref_audio);
    tdcc_scores(j) = max(tdcc);
    
    % fprintf('    Match with %s: %.4f\n', ref_file, tdcc_scores(j)); % Commented out
end

% Normalize scores using self-match (ideal case)
max_dominant_freq_score = max(dominant_freq_scores);
normalized_dominant_freq_scores = dominant_freq_scores / max_dominant_freq_score;

max_spcc_score = max(spcc_scores);
normalized_spcc_scores = spcc_scores / max_spcc_score;

max_tdcc_score = max(tdcc_scores);
normalized_tdcc_scores = tdcc_scores / max_tdcc_score;

% Define weights for each stage (Prioritize TDCC)
dominant_freq_weight = 0.2;  % Lower weight for dominant frequency matching
spcc_weight = 0.2;           % Lower weight for spectrogram cross-correlation
tdcc_weight = 0.6;           % Higher weight for time-domain cross-correlation

% Combine the normalized scores with weights
combined_scores = dominant_freq_weight * normalized_dominant_freq_scores + ...
                  spcc_weight * normalized_spcc_scores + ...
                  tdcc_weight * normalized_tdcc_scores;

% Check if any combined score meets the threshold
[max_combined_score, max_idx] = max(combined_scores);
if max_combined_score >= threshold
    best_match = ref_files{max_idx};
    best_score = combined_scores(max_idx);
    method_used = 'Dominant Frequencies, Spectrogram Cross-Correlation, and Time-Domain Cross-Correlation (Threshold Met)';
    % fprintf('\nBest match determined by %s: %s with score %.4f\n', method_used, best_match, best_score);
    fprintf('\nBest match determined: %s', best_match);

    return;
end

% If no stage meets the threshold, fall back to the best combined score
[~, fallback_idx] = max(combined_scores);
best_match = ref_files{fallback_idx};
best_score = combined_scores(fallback_idx);
method_used = 'Fallback to Best Combined Score (No Threshold Met)';
% fprintf('\nBest match determined by %s: %s with score %.4f\n', method_used, best_match, best_score);
fprintf('\nBest match determined: %s', best_match);


% Helper functions
function [dominant_freqs, dominant_order] = find_dominant_frequencies(spectrogram, frequencies, n)
    [~, idx] = sort(max(spectrogram, [], 2), 'descend');
    dominant_freqs = frequencies(idx(1:n));
    dominant_order = idx(1:n);
end

function similarity_score = compute_similarity_score(ref_dominant_freqs, ref_dominant_order, task_dominant_freqs, task_dominant_order)
    score = 0;
    for i = 1:length(ref_dominant_freqs)
        if ref_dominant_freqs(i) == task_dominant_freqs(i) && ref_dominant_order(i) == task_dominant_order(i)
            score = score + 1;
        end
    end
    similarity_score = score / length(ref_dominant_freqs);
end
