
% Define the audio file and text file
audioFile = '1.wav'; % Replace with your actual audio file path
textFile = '1.txt'; % Replace with your actual text file path

loud_segments = task1(audioFile, textFile);

function [loud_segments] = task1(audioFile, textFile)
    loud_segments = analyzeLoudnessWithNoiseReduction(audioFile, textFile);
end

function [loud_segments] = analyzeLoudnessWithNoiseReduction(audioFile, textFile)
    % Read audio
    [y, fs] = audioread(audioFile);
    cleaned_speech = cleanSpeechAudio(y, fs);
    
    [words, start_times, end_times, loudness_labels] = readTimeFile('1.txt');
    
    % Initialize results
    loud_segments = struct('word', {}, 'start', {}, 'end', {}, 'is_loud', {}, 'energy', {}, 'peak', {}, 'zcr', {});
    
    % Process each segment
    for i = 1:length(words)
        % Convert time to samples
        start_sample = max(1, round(start_times(i) * fs)); % Ensure >= 1
        end_sample = min(length(cleaned_speech), round(end_times(i) * fs)); % Ensure <= signal length
        
        % Validate segment
        if start_sample >= end_sample
            fprintf('Skipping invalid segment: %s [%f - %f]\n', words{i}, start_times(i), end_times(i));
            continue;
        end
        
        % Extract segment from cleaned speech
        segment = cleaned_speech(start_sample:end_sample);
        
        % Calculate features
        rms_energy = rms(segment);
        peak_amp = max(abs(segment));
        zcr = sum(abs(diff(sign(segment)))) / (2 * length(segment));
        
        % Store results in structure
        loud_segments(i).word = words{i};
        loud_segments(i).start = start_times(i);
        loud_segments(i).end = end_times(i);
        loud_segments(i).is_loud = loudness_labels(i);
        loud_segments(i).energy = rms_energy;
        loud_segments(i).peak = peak_amp;
        loud_segments(i).zcr = zcr;
    end
    
    % Visualize results
    plotResults(y, cleaned_speech, fs, loud_segments);
end

function [cleaned_speech] = cleanSpeechAudio(y, fs)
    % Parameters for speech enhancement
    frame_len = round(0.025 * fs);     % 25ms frame
    frame_overlap = round(0.015 * fs);  % 15ms overlap
    
    % 1. Spectral Subtraction for initial noise reduction
    cleaned = spectralSubtraction(y, fs, frame_len, frame_overlap);
    
    % 2. Wiener filtering for further enhancement
    cleaned = wienerFilter(cleaned, fs, frame_len, frame_overlap);
    
    % 3. Apply bandpass filter for speech frequency range (80Hz-3.5kHz)
    cleaned = speechBandpassFilter(cleaned, fs);
    
    cleaned_speech = cleaned;
end

function plotResults(y, cleaned_speech, fs, segments)
    figure('Position', [100 100 1200 800]);
    
    % Plot original waveform
    subplot(2,1,1);
    t = (0:length(y)-1)/fs;
    plot(t, y);
    title('Original Speech Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    ylim([-1.1 1.1]);
    
    % Add word labels and highlight loud segments
    hold on;
    for i = 1:length(segments)
        % Add word label
        text(segments(i).start, 1.0, segments(i).word, ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
        
        % Highlight loud segments
        if segments(i).is_loud == 1
            x = [segments(i).start segments(i).start segments(i).end segments(i).end];
            y = [-1 1 1 -1];
            patch(x, y, 'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    hold off;
    
    % Plot cleaned waveform
    subplot(2,1,2);
    plot(t, cleaned_speech);
    title('Cleaned Speech Signal with Loud Segments Highlighted');
    xlabel('Time (s)');
    ylabel('Amplitude');
    ylim([-1.1 1.1]);
    
    % Add word labels and highlight loud segments
    hold on;
    for i = 1:length(segments)
        % Add word label
        text(segments(i).start, 1.0, segments(i).word, ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
        
        % Highlight loud segments
        if segments(i).is_loud == 1
            x = [segments(i).start segments(i).start segments(i).end segments(i).end];
            y = [-1 1 1 -1];
            patch(x, y, 'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    hold off;
end

% Load audio
[audio, fs] = audioread('1.wav');
[words, start_times, end_times, loudness_indicator] = readTimeFile('1.txt');

for i = 1:length(words)
    fprintf('%s\t\t%f\t%f\t%d\n', words{i}, start_times(i), end_times(i), loudness_indicator(i));
end
fprintf('\n');

% Task 1: With timing
for i = 1:length(words)
    % Convert time to samples
    start_sample = round(start_times(i) * fs);
    end_sample = round(end_times(i) * fs);
    
    % Extract word segment
    word_segment = audio(start_sample:end_sample);
    
    % Calculate energy/loudness features
    energy = rms(word_segment);
    % Add more features as needed
end

% Task 2: Without timing
% Use sliding window
window_size = round(0.02 * fs);  % 20ms window
overlap = round(0.01 * fs);      % 10ms overlap

% Calculate short-time energy
energy = zeros(1, floor((length(audio)-window_size)/overlap));
for i = 1:length(energy)
    start_idx = (i-1)*overlap + 1;
    segment = audio(start_idx:start_idx+window_size-1);
    energy(i) = rms(segment);
end

% Set threshold and detect loud segments
threshold = mean(energy) + std(energy);
loud_segments = energy > threshold;