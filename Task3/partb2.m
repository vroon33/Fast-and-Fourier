function [loudness, time] = computeLoudness(audioFile)
    % Read audio file
    [y, fs] = audioread(audioFile);
    
    % Convert to mono if stereo
    if size(y, 2) > 1
        y = mean(y, 2);
    end
    
    y = noisefilter(y, fs);
    
    % STFT parameters
    n_fft = 2048;
    hop_length = 1024;
    window = hann(n_fft, 'periodic');
    
    % Compute STFT
    [S, f, t] = spectrogram(y, window, n_fft-hop_length, n_fft, fs);
    
    % Convert to magnitude spectrum
    S_mag = abs(S);
    
    % Debug: Print sizes
    disp(['Size of S_mag: ', num2str(size(S_mag))]);
    disp(['Number of frequency points: ', num2str(length(f))]);
    
    % Compute A-weighting coefficients for the actual frequencies
    A = computeAWeighting(f);
    
    % Reshape A-weighting to match spectrogram dimensions
    A = A(:);  % Ensure A is a column vector
    S_weighted = S_mag .* A;  % Broadcasting will handle the multiplication
    
    % Compute RMS for each time frame
    loudness = sqrt(mean(S_weighted.^2, 1));
    
    % Convert to dB scale with reference to max value
    loudness_db = 20 * log10(loudness / max(loudness));
    
    % Plot results
    figure;
    
    % Plot waveform
    subplot(3,1,1);
    plot((0:length(y)-1)/fs, y);
    title('Audio Waveform');
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % Plot spectrogram
    subplot(3,1,2);
    imagesc(t, f, 10*log10(S_mag));
    axis xy;
    title('Spectrogram');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    
    % Plot loudness curve
    subplot(3,1,3);
    plot(t, loudness_db);
    title('A-weighted Loudness');
    xlabel('Time (s)');
    ylabel('Loudness (dB)');
    grid on;
    
    % Add threshold line
    hold on;
    threshold = mean(loudness_db) + 1.5 * std(loudness_db);
    plot([t(1) t(end)], [threshold threshold], 'r--', 'LineWidth', 2);
    
    % Find and mark loud segments
    is_loud = loudness_db > threshold;
    plot(t(is_loud), loudness_db(is_loud), 'g.', 'MarkerSize', 10);
    legend('Loudness', 'Threshold', 'Loud Segments');
    
    % Print loud segments
    transitions = diff([0, is_loud, 0]);
    loud_starts = find(transitions == 1);
    loud_ends = find(transitions == -1) - 1;
    
    fprintf('\nLoud segments detected at:\n');
    for i = 1:length(loud_starts)
        fprintf('%.2f - %.2f seconds\n', t(loud_starts(i)), t(loud_ends(i)));
    end
    
    time = t;
end

function A = computeAWeighting(f)
    % Compute A-weighting coefficients according to IEC 61672:2003
    f = f(:);  % Ensure f is a column vector
    
    % Constants for A-weighting calculation
    f1 = 20.598997; 
    f2 = 107.65265;
    f3 = 737.86223;
    f4 = 12194.217;
    
    % A-weighting formula
    numerator = (f4^2 .* f.^4);
    denominator = ((f2 + f.*1i) .* (f3 + f.*1i) .* (f1 + f.*1i) .* (f4 + f.*1i));
    
    % Compute frequency response
    response = numerator ./ abs(denominator);
    
    % Convert to dB and normalize
    A = 20 * log10(abs(response));
    A = A - max(A); % Normalize
    
    % Convert from dB to linear scale
    A = 10.^(A/20);
end


% Example usage
[loudness, time] = computeLoudness('2.wav');