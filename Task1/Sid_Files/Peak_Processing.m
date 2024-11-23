% MATLAB function to take an audio signal as input and return a vector of 10 highest frequencies and their magnitudes.

function [Top_Peak_Frequencies, Top_Peak_Magnitudes] = Peak_Processing(y, Fs)
    % Peaks atleast 100 Hz apart

    y = noisefilter(y, Fs);

    % Compute FFT
    N = length(y);
    fft_length = N;  % Match the FFT length
    Y = fft(y, fft_length)';
    Y_shifted = fftshift(Y);  % Center the FFT
    mag_spectrum = abs(Y_shifted);  % Get magnitude spectrum

    % Create frequency vector centered at 0
    freq = (-fft_length/2:fft_length/2-1) * (Fs/fft_length);

    % Find peaks with a minimum distance of 100 Hz
    min_distance = 100 / (Fs / length(freq));  % Convert 100 Hz to index distance
    [peaks, peak_locs] = findpeaks(mag_spectrum, 'MinPeakDistance', min_distance);

    % Sort peaks by magnitude (descending order) and get the top 5
    [sorted_peaks, sorted_indices] = sort(peaks, 'descend');
    top_peaks = sorted_peaks(1:min(20, length(sorted_peaks)));  % Take top 10 peaks
    top_peak_locs = peak_locs(sorted_indices(1:min(20, length(sorted_peaks))));
    top_peak_freqs = freq(top_peak_locs);  % Corresponding frequencies

    Top_Peak_Frequencies = top_peak_freqs(2:2:end);
    Top_Peak_Magnitudes = top_peaks(2:2:end);

end