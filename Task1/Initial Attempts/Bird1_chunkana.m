% Function to plot FFT of a specific chunk
function plotChunkFFT(chunk_idx, chunk_samples, overlap, y, Fs)
    % Extract the specified chunk
    start_idx = (chunk_idx-1)*(chunk_samples-overlap) + 1;
    end_idx = min(start_idx + chunk_samples - 1, length(y));
    chunk = y(start_idx:end_idx);
    
    % Pad chunk to power of 2
    chunk = [chunk; zeros(2^nextpow2(length(chunk)) - length(chunk), 1)];
    N = length(chunk);
    
    % Compute FFT
    Y = fft(chunk);
    Y_shifted = fftshift(Y);
    mag_spectrum = abs(Y_shifted);
    
    % Create frequency vector
    freq = (-N/2:N/2-1)*(Fs/N);
    
    % Consider positive frequencies
    positive_freq_idx = freq > 0;
    positive_mag_spectrum = mag_spectrum(positive_freq_idx);
    positive_freq = freq(positive_freq_idx);
    
    % Find peaks
    [peaks, locs] = findpeaks(positive_mag_spectrum, 'MinPeakProminence', max(positive_mag_spectrum)/20);
    peak_freqs = positive_freq(locs);
    
    % Create the plot
    figure;
    plot(positive_freq, positive_mag_spectrum, 'b-');
    hold on;
    plot(peak_freqs, peaks, 'ro', 'MarkerSize', 10);
    
    % Add labels and title
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(sprintf('FFT of Chunk %d (%.2f-%.2f seconds)', ...
        chunk_idx, ...
        (chunk_idx-1)*chunk_samples/Fs, ...
        min(chunk_idx*chunk_samples/Fs, length(y)/Fs)));
    grid on;
    
    % Add text annotations for peak frequencies
    for i = 1:min(5, length(peak_freqs))
        text(peak_freqs(i), peaks(i), ...
            sprintf('%.1f Hz', peak_freqs(i)), ...
            'VerticalAlignment', 'bottom');
    end
end


% Read the audio file
[y, Fs] = audioread('Bird1.wav');

% Define chunk parameters
chunk_duration = 0.1; % seconds
chunk_samples = round(chunk_duration * Fs);
overlap = 0; % Can be adjusted if needed

% Calculate number of chunks
num_chunks = floor(length(y)/(chunk_samples - overlap));

% Initialize storage for results
all_dominant_freqs = cell(num_chunks, 1);

% Process each chunk
for chunk_idx = 1:num_chunks
    % Extract chunk
    start_idx = (chunk_idx-1)*(chunk_samples-overlap) + 1;
    end_idx = min(start_idx + chunk_samples - 1, length(y));
    chunk = y(start_idx:end_idx);
    
    % Pad chunk to power of 2
    chunk = [chunk; zeros(2^nextpow2(length(chunk)) - length(chunk), 1)];
    N = length(chunk);
    
    % Compute FFT
    Y = fft(chunk);
    Y_shifted = fftshift(Y);
    mag_spectrum = abs(Y_shifted);
    
    % Create frequency vector
    freq = (-N/2:N/2-1)*(Fs/N);
    
    % Consider positive frequencies
    positive_freq_idx = freq > 0;
    positive_mag_spectrum = mag_spectrum(positive_freq_idx);
    positive_freq = freq(positive_freq_idx);
    
    % Find peaks
    [peaks, locs] = findpeaks(positive_mag_spectrum, 'MinPeakProminence', max(positive_mag_spectrum)/20);
    peak_freqs = positive_freq(locs);
    
    % Sort peaks
    [sorted_peaks, idx] = sort(peaks, 'descend');
    dominant_freqs = peak_freqs(idx);
    
    % Store first 5 dominant frequencies for this chunk
    all_dominant_freqs{chunk_idx} = dominant_freqs(1:min(5, length(dominant_freqs)));
    
    % Optional: Display time and frequencies for each chunk
    time_start = (chunk_idx-1)*chunk_duration;
    fprintf('Chunk %d (%.2f-%.2f s): Dominant frequencies = %.1f Hz, %.1f Hz, %.1f Hz\n', ...
        chunk_idx, time_start, time_start+chunk_duration, ...
        all_dominant_freqs{chunk_idx}(1:min(3, length(all_dominant_freqs{chunk_idx}))));
end
plotChunkFFT(1, chunk_samples, overlap, y, Fs);


% Example usage (add this to your main code):
% plotChunkFFT(1, chunk_samples, overlap, y, Fs);  % Plot first chunk