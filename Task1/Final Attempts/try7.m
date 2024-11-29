% List the reference audio files
ref_files = {'bird1.wav', 'bird2.wav', 'bird3.wav'};

% Load the task audio file
task_file = 'F8.wav';
[task_audio, task_fs] = audioread(task_file);

% Initialize variables to track scores at each stage
dominant_freq_scores = zeros(1, length(ref_files));
spcc_scores = zeros(1, length(ref_files));
tdcc_scores = zeros(1, length(ref_files));

% Iterate through the reference audio files
for j = 1:length(ref_files)
    ref_file = ref_files{j};
    
    % Load the reference audio file
    [ref_audio, ref_fs] = audioread(ref_file);
    
    % Compute the spectrograms
    [task_spectrogram, task_frequencies, ~] = spectrogram(task_audio, hamming(256), 128, 256, task_fs);
    [ref_spectrogram, ref_frequencies, ~] = spectrogram(ref_audio, hamming(256), 128, 256, ref_fs);
    
    % Stage 1: Compute scores based on dominant frequencies
    [task_dominant_freqs, task_dominant_order] = find_dominant_frequencies(task_spectrogram, task_frequencies, 4);
    [ref_dominant_freqs, ref_dominant_order] = find_dominant_frequencies(ref_spectrogram, ref_frequencies, 4);
    dominant_freq_scores(j) = compute_similarity_score(ref_dominant_freqs, ref_dominant_order, task_dominant_freqs, task_dominant_order);
end

% Analyze dominant frequency scores
[max_dom_score, max_idx] = max(dominant_freq_scores);
sorted_dom_scores = sort(dominant_freq_scores, 'descend');
if length(sorted_dom_scores) > 1
    score_difference = max_dom_score - sorted_dom_scores(2);
else
    score_difference = max_dom_score; % Single score case
end
variance_dom = var(dominant_freq_scores);

% Check if the dominant frequency stage is decisive
if score_difference > 0.2 || variance_dom > 0.1
    best_match = ref_files{max_idx};
    fprintf('Best match based on dominant frequencies: %s\n', best_match);
    return;
end

% Stage 2: Spectrogram cross-correlation
for j = 1:length(ref_files)
    [ref_audio, ref_fs] = audioread(ref_files{j});
    [ref_spectrogram, ~, ~] = spectrogram(ref_audio, hamming(256), 128, 256, ref_fs);
    [spcc, ~] = xcorr(task_spectrogram(:), ref_spectrogram(:));
    spcc_scores(j) = max(spcc);
end

% Analyze spectrogram cross-correlation scores
[max_spcc_score, max_idx] = max(spcc_scores);
sorted_spcc_scores = sort(spcc_scores, 'descend');
if length(sorted_spcc_scores) > 1
    score_difference = max_spcc_score - sorted_spcc_scores(2);
else
    score_difference = max_spcc_score; % Single score case
end
variance_spcc = var(spcc_scores);

% Check if the spectrogram cross-correlation stage is decisive
if score_difference > 0.15 || variance_spcc > 0.1
    best_match = ref_files{max_idx};
    fprintf('Best match based on spectrogram cross-correlation: %s\n', best_match);
    return;
end

% Stage 3: Time-domain cross-correlation
for j = 1:length(ref_files)
    [ref_audio, ~] = audioread(ref_files{j});
    [tdcc, ~] = xcorr(task_audio, ref_audio);
    tdcc_scores(j) = max(tdcc);
end

% Analyze time-domain cross-correlation scores
[max_tdcc_score, max_idx] = max(tdcc_scores);
sorted_tdcc_scores = sort(tdcc_scores, 'descend');
if length(sorted_tdcc_scores) > 1
    score_difference = max_tdcc_score - sorted_tdcc_scores(2);
else
    score_difference = max_tdcc_score; % Single score case
end
variance_tdcc = var(tdcc_scores);

% Check if the time-domain cross-correlation stage is decisive
if score_difference > 0.1 || variance_tdcc > 0.05
    best_match = ref_files{max_idx};
    fprintf('Best match based on time-domain cross-correlation: %s\n', best_match);
    return;
end

% If no stage is decisive, output the highest overall match
fprintf('No clear match found; best match overall is %s\n', ref_files{max_idx});

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
