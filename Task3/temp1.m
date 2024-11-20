% task1.m
function [loud_segments] = task1(audioFile, textFile)
    loud_segments = analyzeLoudnessWithNoiseReduction(audioFile, textFile);
end

function [loud_segments] = analyzeLoudnessWithNoiseReduction(audioFile, textFile)
    % Read audio
    [y, fs] = audioread(audioFile);
    cleaned_speech = cleanSpeechAudio(y, fs);
    
    [words, start_times, end_times, loudness_indicator] = readTimeFile(textFile);
    
    % Initialize results
    loud_segments = struct('word', {}, 'start', {}, 'end', {}, 'is_loud', {}, 'energy', {}, 'peak', {}, 'zcr', {});
    
    % Process each segment
    for i = 1:length(words)
        % Convert time to samples
        start_sample = max(1, round(start_times(i) * fs));
        end_sample = min(length(cleaned_speech), round(end_times(i) * fs));
        
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
        loud_segments(i).is_loud = loudness_indicator(i);
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
    
    % 1. Spectral Subtraction
    cleaned = spectralSubtraction(y, fs, frame_len, frame_overlap);
    
    % 2. Wiener filtering
    cleaned = wienerFilter(cleaned, fs, frame_len, frame_overlap);
    
    % 3. Apply bandpass filter for speech frequency range (80Hz-3.5kHz)
    cleaned = speechBandpassFilter(cleaned, fs);
    
    cleaned_speech = cleaned;
end

function [enhanced_signal] = spectralSubtraction(y, fs, frame_len, frame_overlap)
    % Implementation of spectral subtraction
    hop_len = frame_len - frame_overlap;
    num_frames = floor((length(y) - frame_overlap) / hop_len);
    window = hamming(frame_len);
    
    % Estimate noise from first few frames (assumed to be noise-only)
    noise_frames = 5;
    noise_spectrum = zeros(frame_len, 1);
    
    for i = 1:noise_frames
        frame = y((i-1)*hop_len + 1 : (i-1)*hop_len + frame_len);
        noise_spectrum = noise_spectrum + abs(fft(frame .* window));
    end
    noise_spectrum = noise_spectrum / noise_frames;
    
    % Process each frame
    enhanced = zeros(size(y));
    for i = 1:num_frames
        % Extract frame
        start_idx = (i-1)*hop_len + 1;
        frame = y(start_idx : min(start_idx + frame_len - 1, length(y)));
        
        if length(frame) < frame_len
            frame(end+1:frame_len) = 0;
        end
        
        % Apply window
        windowed_frame = frame .* window;
        
        % FFT
        spec = fft(windowed_frame);
        mag_spec = abs(spec);
        phase_spec = angle(spec);
        
        % Spectral subtraction
        subtracted_spec = max(mag_spec - noise_spectrum, 0.01 * mag_spec);
        
        % Reconstruct signal
        enhanced_frame = real(ifft(subtracted_spec .* exp(1j * phase_spec)));
        
        % Overlap-add
        enhanced(start_idx : start_idx + frame_len - 1) = ...
            enhanced(start_idx : start_idx + frame_len - 1) + ...
            enhanced_frame(1:min(frame_len, length(y) - start_idx + 1));
    end
    
    enhanced_signal = enhanced;
end

function [enhanced_signal] = wienerFilter(y, fs, frame_len, frame_overlap)
    % Implementation of Wiener filter
    hop_len = frame_len - frame_overlap;
    num_frames = floor((length(y) - frame_overlap) / hop_len);
    window = hamming(frame_len);
    
    % Estimate noise from first few frames
    noise_frames = 5;
    noise_psd = zeros(frame_len, 1);
    
    for i = 1:noise_frames
        frame = y((i-1)*hop_len + 1 : (i-1)*hop_len + frame_len);
        noise_psd = noise_psd + abs(fft(frame .* window)).^2;
    end
    noise_psd = noise_psd / noise_frames;
    
    % Process each frame
    enhanced = zeros(size(y));
    for i = 1:num_frames
        % Extract frame
        start_idx = (i-1)*hop_len + 1;
        frame = y(start_idx : min(start_idx + frame_len - 1, length(y)));
        
        if length(frame) < frame_len
            frame(end+1:frame_len) = 0;
        end
        
        % Apply window
        windowed_frame = frame .* window;
        
        % FFT
        spec = fft(windowed_frame);
        sig_psd = abs(spec).^2;
        phase_spec = angle(spec);
        
        % Wiener filter
        H = max(1 - noise_psd ./ (sig_psd + eps), 0.1);
        enhanced_spec = abs(spec) .* H;
        
        % Reconstruct signal
        enhanced_frame = real(ifft(enhanced_spec .* exp(1j * phase_spec)));
        
        % Overlap-add
        enhanced(start_idx : start_idx + frame_len - 1) = ...
            enhanced(start_idx : start_idx + frame_len - 1) + ...
            enhanced_frame(1:min(frame_len, length(y) - start_idx + 1));
    end
    
    enhanced_signal = enhanced;
end

function [filtered_signal] = speechBandpassFilter(y, fs)
    % Design bandpass filter for speech (80Hz-3.5kHz)
    low_freq = 80;
    high_freq = 3500;
    
    [b, a] = butter(4, [low_freq high_freq]/(fs/2), 'bandpass');
    filtered_signal = filtfilt(b, a, y);
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
    
    % Add word labels, highlight loud segments, and vertical lines
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
        
        % Add vertical lines at word boundaries
        plot([segments(i).start segments(i).start], [-1 1], ':', 'Color', [0.5 0.5 0.5]);
        plot([segments(i).end segments(i).end], [-1 1], ':', 'Color', [0.5 0.5 0.5]);
    end
    hold off;
    
    % Plot cleaned waveform
    subplot(2,1,2);
    plot(t, cleaned_speech);
    title('Cleaned Speech Signal with Loud Segments Highlighted');
    xlabel('Time (s)');
    ylabel('Amplitude');
    ylim([-1.1 1.1]);
    
    % Add word labels, highlight loud segments, and vertical lines
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
        
        % Add vertical lines at word boundaries
        plot([segments(i).start segments(i).start], [-1 1], ':', 'Color', [0.5 0.5 0.5]);
        plot([segments(i).end segments(i).end], [-1 1], ':', 'Color', [0.5 0.5 0.5]);
    end
    hold off;
end

function [words, start_times, end_times, loudness_indicator] = readTimeFile(filename)
    % Read the text file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Initialize arrays
    words = {};
    start_times = [];
    end_times = [];
    loudness_indicator = [];
    
    % Read line by line
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line)
            % Parse line (assuming format: word start_time end_time is_loud)
            parts = strsplit(line);
            words{end+1} = parts{1};
            start_times(end+1) = str2double(parts{2});
            end_times(end+1) = str2double(parts{3});
            loudness_indicator(end+1) = str2double(parts{4});
        end
    end
    
    fclose(fid);
end