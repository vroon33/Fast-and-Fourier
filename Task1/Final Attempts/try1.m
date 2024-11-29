% List the reference audio files
ref_files = {'bird1.wav', 'bird2.wav', 'bird3.wav'};

% Load the task audio file
task_file = 'F8.wav';
[task_audio, task_fs] = audioread(task_file);

% Initialize variables to track the best match
best_match = '';
best_correlation = -inf;

% Iterate through the reference audio files
for j = 1:length(ref_files)
    ref_file = ref_files{j};
    
    % Load the reference audio file
    [ref_audio, ref_fs] = audioread(ref_file);
    
    % Compute the spectrograms
    [task_spectrogram, task_frequencies, task_time] = spectrogram(task_audio, hamming(256), 128, 256, task_fs);
    [ref_spectrogram, ref_frequencies, ref_time] = spectrogram(ref_audio, hamming(256), 128, 256, ref_fs);
    
    % Compute the spectrogram cross-correlation (SPCC)
    [spcc, lags] = xcorr(task_spectrogram(:), ref_spectrogram(:));
    
    % Find the maximum correlation and its index
    [max_corr, max_idx] = max(spcc);
    
    % Update the best match if the current correlation is higher
    if max_corr > best_correlation
        best_match = ref_file;
        best_correlation = max_corr;
    end
end

% Display the best match for the task audio file
fprintf('Best match for %s is %s\n', task_file, best_match);