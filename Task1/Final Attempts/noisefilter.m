function [cleaned_signal] = noisefilter(input_signal, fs)
    % Convert stereo to mono if necessary
    if size(input_signal, 2) > 1
        input_signal = mean(input_signal, 2);
    end
    
    % Ensure input_signal is a column vector
    input_signal = input_signal(:);
    
    % Parameters
    frame_length = round(0.025 * fs);    % 25ms frame
    frame_shift = round(0.010 * fs);     % 10ms shift
    fft_length = 2^nextpow2(frame_length);
    
    % Initialize noise estimation from first few frames (assumed to be noise)
    num_noise_frames = 5;
    noise_estimate = estimate_noise(input_signal(1:frame_length*num_noise_frames), ...
        frame_length, frame_shift, fft_length);
    
    % Frame the signal
    num_frames = floor((length(input_signal) - frame_length) / frame_shift) + 1;
    cleaned_frames = zeros(frame_length, num_frames);
    
    % Create window once
    window = hanning(frame_length);
    
    % Process each frame
    for i = 1:num_frames
        % Extract frame
        start_idx = (i-1) * frame_shift + 1;
        frame = input_signal(start_idx:start_idx + frame_length - 1);
        
        % Apply Hanning window
        windowed_frame = frame .* window;
        
        % FFT
        spectrum = fft(windowed_frame, fft_length);
        mag_spectrum = abs(spectrum);
        phase_spectrum = angle(spectrum);
        
        % Spectral subtraction
        alpha = 2.0;    % Over-subtraction factor
        beta = 0.01;    % Spectral floor
        
        % Compute enhanced magnitude spectrum
        enhanced_mag = max(mag_spectrum - alpha * noise_estimate, ...
            beta * mag_spectrum);
        
        % Reconstruct signal
        enhanced_spectrum = enhanced_mag .* exp(1j * phase_spectrum);
        enhanced_frame = real(ifft(enhanced_spectrum, fft_length));
        enhanced_frame = enhanced_frame(1:frame_length);
        
        % Ensure enhanced_frame is a column vector
        enhanced_frame = enhanced_frame(:);
        
        % Save frame
        cleaned_frames(:,i) = enhanced_frame .* window;
    end
    
    % Overlap-add synthesis
    cleaned_signal = zeros(length(input_signal), 1);
    for i = 1:num_frames
        start_idx = (i-1) * frame_shift + 1;
        cleaned_signal(start_idx:start_idx + frame_length - 1) = ...
            cleaned_signal(start_idx:start_idx + frame_length - 1) + ...
            cleaned_frames(:,i);
    end
    
    % Normalize output
    cleaned_signal = cleaned_signal / max(abs(cleaned_signal));
end

function noise_estimate = estimate_noise(noise_segment, frame_length, frame_shift, fft_length)
    % Ensure noise_segment is a column vector
    noise_segment = noise_segment(:);
    
    % Frame the noise segment
    num_frames = floor((length(noise_segment) - frame_length) / frame_shift) + 1;
    noise_power = zeros(fft_length, 1);
    
    % Create window once
    window = hanning(frame_length);
    
    % Process each frame
    for i = 1:num_frames
        start_idx = (i-1) * frame_shift + 1;
        frame = noise_segment(start_idx:start_idx + frame_length - 1);
        
        % Window and FFT
        windowed_frame = frame .* window;
        spectrum = fft(windowed_frame, fft_length);
        noise_power = noise_power + abs(spectrum);
    end
    
    % Average noise power spectrum
    noise_estimate = noise_power / num_frames;
end