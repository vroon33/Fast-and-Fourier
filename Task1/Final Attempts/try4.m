% List the reference audio files
ref_files = {'bird1.wav', 'bird2.wav', 'bird3.wav'};

% Load the task audio file
task_file = 'F8.wav';
[task_audio, task_fs] = audioread(task_file);

% Initialize variables to track the best match
best_match = '';
best_score = 0;

% Iterate through the reference audio files
for j = 1:length(ref_files)
    ref_file = ref_files{j};
    
    % Load the reference audio file
    [ref_audio, ref_fs] = audioread(ref_file);
    
    % Compute the spectrograms
    [task_spectrogram, task_frequencies, task_time] = spectrogram(task_audio, hamming(256), 128, 256, task_fs);
    [ref_spectrogram, ref_frequencies, ref_time] = spectrogram(ref_audio, hamming(256), 128, 256, ref_fs);
    
    % Find the 4 dominant frequencies and their order in the spectrograms
    [task_dominant_freqs, task_dominant_order] = find_dominant_frequencies(task_spectrogram, task_frequencies, 4);
    [ref_dominant_freqs, ref_dominant_order] = find_dominant_frequencies(ref_spectrogram, ref_frequencies, 4);
    
    % Compute a similarity score based on the dominant frequencies and their order
    similarity_score = compute_similarity_score(ref_dominant_freqs, ref_dominant_order, task_dominant_freqs, task_dominant_order);
    
    % Update the best match if the current score is higher
    if similarity_score > best_score
        best_match = ref_file;
        best_score = similarity_score;
    end
end

% Display the best match for the task audio file
fprintf('Best match for %s is %s\n', task_file, best_match);

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