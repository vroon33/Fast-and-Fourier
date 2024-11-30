% Read and analyze reference files
reference_folder = 'Reference';
ref_files = dir(fullfile(reference_folder, '*.wav'));
num_ref = length(ref_files);

% Initialize matrices to store features
ref_features = zeros(num_ref, 4);

% Extract features from reference files
for i = 1:num_ref
    [y, fs] = audioread(fullfile(reference_folder, ref_files(i).name));
    
    % Convert to mono if stereo
    if size(y, 2) > 1
        y = mean(y, 2);
    end
    
    % Calculate features
    % 1. Average magnitude
    ref_features(i, 1) = mean(abs(y));
    
    % 2. Spectral centroid
    window = hamming(round(fs*0.03));
    noverlap = round(length(window)*0.5);
    [~, f, t, p] = spectrogram(y, window, noverlap, [], fs);
    centroid = sum(f .* p) ./ sum(p);
    ref_features(i, 2) = mean(centroid);
    
    % 3. Zero crossing rate
    ref_features(i, 3) = sum(abs(diff(sign(y)))) / (2*length(y));
    
    % 4. Spectral rolloff
    cumsum_spec = cumsum(p);
    rolloff_point = 0.85 * sum(p);
    [~, rolloff_idx] = min(abs(cumsum_spec - rolloff_point));
    ref_features(i, 4) = f(rolloff_idx);
end

% Process task files
task_folder = 'Task';
task_files = dir(fullfile(task_folder, '*.wav'));
num_task = length(task_files);

% Initialize matrix for task features
task_features = zeros(num_task, 4);

% Extract features from task files
for i = 1:num_task
    [y, fs] = audioread(fullfile(task_folder, task_files(i).name));
    
    % Convert to mono if stereo
    if size(y, 2) > 1
        y = mean(y, 2);
    end
    
    % Calculate same features as reference
    % 1. Average magnitude
    task_features(i, 1) = mean(abs(y));
    
    % 2. Spectral centroid
    window = hamming(round(fs*0.03));
    noverlap = round(length(window)*0.5);
    [~, f, t, p] = spectrogram(y, window, noverlap, [], fs);
    centroid = sum(f .* p) ./ sum(p);
    task_features(i, 2) = mean(centroid);
    
    % 3. Zero crossing rate
    task_features(i, 3) = sum(abs(diff(sign(y)))) / (2*length(y));
    
    % 4. Spectral rolloff
    cumsum_spec = cumsum(p);
    rolloff_point = 0.85 * sum(p);
    [~, rolloff_idx] = min(abs(cumsum_spec - rolloff_point));
    task_features(i, 4) = f(rolloff_idx);
end

% Classify each task file using correlation with reference features
predictions = zeros(num_task, 1);
for i = 1:num_task
    corr_scores = zeros(num_ref, 1);
    for j = 1:num_ref
        % Calculate correlation between feature vectors
        corr_matrix = corrcoef(task_features(i,:), ref_features(j,:));
        corr_scores(j) = corr_matrix(1,2);
    end
    % Assign class based on highest correlation
    [~, predictions(i)] = max(corr_scores);
end

% Display results
for i = 1:num_task
    fprintf('Task file %d matches Reference bird species %d\n', i, predictions(i));
end