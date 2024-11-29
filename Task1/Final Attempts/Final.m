ref_files = {'bird1.wav', 'bird2.wav', 'bird3.wav'};

task_file = 'F8.wav';
[task_audio, task_fs] = audioread(task_file);

best_match = '';
best_score = -inf;
method_used = '';

% Initialize score arrays
dominant_freq_scores = zeros(1, length(ref_files));
spcc_scores = zeros(1, length(ref_files));
tdcc_scores = zeros(1, length(ref_files));

% Stage 1: Dominant Frequency Matching + Spectrogram Cross-Correlation (SPCC)
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
end

% Stage 2: Time-Domain Cross-Correlation (TDCC)
for j = 1:length(ref_files)
    ref_file = ref_files{j};
    [ref_audio, ~] = audioread(ref_file);
    
    % Compute TDCC
    [tdcc, ~] = xcorr(task_audio, ref_audio);
    tdcc_scores(j) = max(tdcc);
end

% Normalize scores using self-match (ideal case)
% Ensure that max score is not zero before normalizing
max_dominant_freq_score = max(dominant_freq_scores);
if max_dominant_freq_score > 0
    normalized_dominant_freq_scores = dominant_freq_scores / max_dominant_freq_score;
else
    normalized_dominant_freq_scores = dominant_freq_scores; % No normalization if max score is 0
end

max_spcc_score = max(spcc_scores);
if max_spcc_score > 0
    normalized_spcc_scores = spcc_scores / max_spcc_score;
else
    normalized_spcc_scores = spcc_scores; % No normalization if max score is 0
end

max_tdcc_score = max(tdcc_scores);
if max_tdcc_score > 0
    normalized_tdcc_scores = tdcc_scores / max_tdcc_score;
else
    normalized_tdcc_scores = tdcc_scores; % No normalization if max score is 0
end

% Define weights for each stage (Prioritize TDCC)
dominant_freq_weight = 0.25;  % Lower weight for dominant frequency matching
spcc_weight = 0.35;           % Lower weight for spectrogram cross-correlation
tdcc_weight = 0.4;           % Higher weight for time-domain cross-correlation

% Combine the normalized scores with weights
combined_scores = dominant_freq_weight * normalized_dominant_freq_scores + ...
                  spcc_weight * normalized_spcc_scores + ...
                  tdcc_weight * normalized_tdcc_scores;

% Fall back to the best combined score
[~, best_idx] = max(combined_scores);
best_match = ref_files{best_idx};
best_score = combined_scores(best_idx);
method_used = 'Fallback to Best Combined Score (No Threshold)';

fprintf('\nBest match determined by %s: %s with score %.4f\n', method_used, best_match, best_score);

% Load the best match and compute its spectrogram
[best_match_audio, best_match_fs] = audioread(best_match);
[task_spectrogram, task_frequencies, ~] = spectrogram(task_audio, hamming(256), 128, 256, task_fs);
[best_match_spectrogram, best_match_frequencies, ~] = spectrogram(best_match_audio, hamming(256), 128, 256, best_match_fs);

% Plot the spectrograms of the task file and best match using subplot
figure;
subplot(1, 2, 1);
imagesc(task_frequencies, 1:size(task_spectrogram, 2), abs(task_spectrogram));
axis xy;
title('Task Audio Spectrogram');
xlabel('Frequency (Hz)');
ylabel('Time (s)');
colorbar;

subplot(1, 2, 2);
imagesc(best_match_frequencies, 1:size(best_match_spectrogram, 2), abs(best_match_spectrogram));
axis xy;
title('Best Match Spectrogram');
xlabel('Frequency (Hz)');
ylabel('Time (s)');
colorbar;

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
