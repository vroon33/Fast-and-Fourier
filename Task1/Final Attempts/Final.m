% Bird Matcher

clc, clearvars;

ref_files = {'bird1.wav', 'bird2.wav', 'bird3.wav'};

task_file = 'F6.wav'; % Change the Test file here

[task_audio, task_fs] = audioread(task_file);
% task_audio = flip(task_audio);
% sound(task_audio, task_fs);

% Initialize score arrays
dominant_freq_scores = zeros(1, length(ref_files));
spcc_scores = zeros(1, length(ref_files)); % Spectrogram Cross-Correlation
tdcc_scores = zeros(1, length(ref_files)); % Time Domain Cross-Correlation

% Dominant Frequency Matching + Spectrogram Cross-Correlation (SPCC)
for j = 1:length(ref_files)
    ref_file = ref_files{j};
    [ref_audio, ref_fs] = audioread(ref_file);
    
    % Compute spectrograms, using a Hamming Window for the STFTs
    [task_spectrogram, task_frequencies, ~] = spectrogram(task_audio, hamming(256), 128, 256, task_fs, 'yaxis');
    [ref_spectrogram, ref_frequencies, ~] = spectrogram(ref_audio, hamming(256), 128, 256, ref_fs, 'yaxis');
    
    % Compute dominant frequencies and their order
    [task_dominant_freqs, task_dominant_order] = find_dominant_frequencies(task_spectrogram, task_frequencies, 4);
    [ref_dominant_freqs, ref_dominant_order] = find_dominant_frequencies(ref_spectrogram, ref_frequencies, 4);
    
    % Compute similarity score based on dominant frequencies
    dominant_freq_scores(j) = compute_similarity_score(ref_dominant_freqs, ref_dominant_order, task_dominant_freqs, task_dominant_order);
    
    % Compute spectrogram cross-correlation
    [spcc, ~] = xcorr(task_spectrogram(:), ref_spectrogram(:));
    spcc_scores(j) = max(spcc);
end

% Time-Domain Cross-Correlation (TDCC)
for j = 1:length(ref_files)
    ref_file = ref_files{j};
    [ref_audio, ~] = audioread(ref_file);
    
    % Compute TDCC
    [tdcc, ~] = xcorr(task_audio, ref_audio);
    tdcc_scores(j) = max(tdcc);
end

max_dominant_freq_score = 4; % The max score that a file can have is 4, when all 4 dominant frequencies match
normalized_dominant_freq_scores = dominant_freq_scores / max_dominant_freq_score;

% Normalize scores using self-match (ideal case)
normalized_spcc_scores = zeros(1, length(spcc_scores));
normalized_tdcc_scores = zeros(1, length(tdcc_scores));
max_spcc_score = zeros(1, length(normalized_tdcc_scores));
max_tdcc_score = zeros(1, length(normalized_tdcc_scores));

for ref_bird = 1 : length(ref_files)
    % Ensure that max score is not zero before normalizing
    
    [audio_Data, audio_fs] = audioread(ref_files{ref_bird});
    bird_spect = spectrogram(audio_Data, hamming(256), 128, 256, audio_fs, 'yaxis'); 
    max_spcc_score(ref_bird) = max(xcorr(bird_spect(:) , bird_spect(:)));
    if max_spcc_score(ref_bird) > 0
        normalized_spcc_scores(ref_bird) = spcc_scores(ref_bird) / max_spcc_score(ref_bird);
    else
        normalized_spcc_scores(ref_bird) = spcc_scores(ref_bird); % No normalization if max score is 0
    end

    max_tdcc_score(ref_bird) = max(xcorr(audio_Data , audio_Data));
    if max_tdcc_score(ref_bird) > 0
        normalized_tdcc_scores(ref_bird) = tdcc_scores(ref_bird) / max_tdcc_score(ref_bird);
    else
        normalized_tdcc_scores(ref_bird) = tdcc_scores(ref_bird); % No normalization if max score is 0
    end
end

% Define weights for each stage
dominant_freq_weight = 0.25;  % Lower weight for dominant frequency matching
spcc_weight = 0.35;           % Lower weight for spctrogram cross-correlation
tdcc_weight = 0.4;           % Higher weight for time-domain cross-correlation

% Combine the normalized scores with weights
combined_scores = dominant_freq_weight * normalized_dominant_freq_scores ...
                  + spcc_weight * normalized_spcc_scores ...
                  + tdcc_weight * normalized_tdcc_scores;

% Choose the best combined score
[best_score, best_idx] = max(combined_scores);
best_match = ref_files{best_idx};
[best_match_audio, best_match_fs] = audioread(best_match);

% fprintf('\nBird 1 Match : %.4f', combined_scores(1));
% fprintf('\nBird 2 Match : %.4f', combined_scores(2));
% fprintf('\nBird 3 Match : %.4f\n', combined_scores(3));
fprintf('\nBest match : %s ', best_match);

% Plot the spectrograms of the task file and best match using subplot
figure;

% Plot the best match spectrogram
subplot(2, 1, 2);
spectrogram(best_match_audio, hamming(256), 128, 256, best_match_fs, 'yaxis');
axis xy;
colormap('jet');
title(['Best Match Spectrogram (' best_match ')']);
xlabel('Time (s)');
ylabel('Frequency (kHz)');
colorbar;

% Plot the task audio spectrogram
subplot(2, 1, 1);
spectrogram(task_audio, hamming(256), 128, 256, task_fs, 'yaxis');
axis xy;
colormap('jet');
title(['Task Audio Spectrogram (' task_file ')']);
xlabel('Time (s)');
ylabel('Frequency (kHz)');
colorbar;

function [dominant_freqs, dominant_order] = find_dominant_frequencies(spectrogram, frequencies, n)
    % Returns a vector of the dominant frequencies and their indices in the frequencies array
    [~, idx] = sort(max(spectrogram, [], 2), 'descend'); % Gives us the indices of the max frequency components
    dominant_freqs = frequencies(idx(1:n));
    dominant_order = idx(1:n);
end

function similarity_score = compute_similarity_score(ref_dominant_freqs, ref_dominant_order, task_dominant_freqs, task_dominant_order)
    score = 0;
    for i = 1:length(ref_dominant_freqs)
        if abs(ref_dominant_freqs(i) - task_dominant_freqs(i)) < 100
            score = score + 1;
        end
    end
    similarity_score = score / length(ref_dominant_freqs);
end