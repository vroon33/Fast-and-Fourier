<!-- This is document for organising the codes separated randomly across files -->
# Task 3 (Loudness Detection)

## Parameters used for Loudness Analysis
- Peak Amplitude
- Energy
- Normalized Energy
- Band Energy (energy in fixed frequency bands)
- Normalized Band Energy

## Functions Utilized
- **Noise Filter**

 ```matlab
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
    noise_estimate = estimate_noise(input_signal(1:frame_length*num_noise_frames), frame_length, frame_shift, fft_length);
    
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
        enhanced_mag = max(mag_spectrum - alpha * noise_estimate, beta * mag_spectrum);
        
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
        cleaned_signal(start_idx:start_idx + frame_length - 1) = cleaned_signal(start_idx:start_idx + frame_length - 1) + cleaned_frames(:,i);
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
 ```

#### Code explanation:

The noise filter implementation consists of two main functions: `noisefilter` and `estimate_noise`. Here's a detailed line-by-line explanation:

Main Function (noisefilter):
1. Input Processing [Lines 2-8]:
   - Line 3-5: Converts stereo to mono by averaging channels if necessary
   - Lines 7-8: Ensures input signal is in column vector format

2. Parameters Setup [Lines 10-13]:
   - Line 11: Frame length set to 25ms windows of the signal
   - Line 12: Frame shift set to 10ms overlap between frames
   - Line 13: FFT length calculated as power of 2 for efficient computation

3. Noise Estimation [Lines 15-17]:
   - Line 16: Sets first 5 frames for noise estimation
   - Line 17: Calls estimate_noise function for initial noise profile

4. Frame Processing [Lines 19-33]:
   - Lines 20-21: Segments signal into overlapping frames
   - Line 24: Creates Hanning window for spectral analysis
   - Lines 28-33: Applies windowing and performs FFT
   - Lines 32-33: Extracts magnitude and phase spectra

5. Spectral Subtraction [Lines 35-40]:
   - Lines 36-37: Defines oversubtraction (α=2.0) and floor (β=0.01) parameters
   - Line 40: Implements spectral subtraction formula

6. Signal Reconstruction [Lines 42-59]:
   - Lines 42-45: Reconstructs enhanced frames using IFFT
   - Lines 54-58: Performs overlap-add synthesis
   - Line 59: Normalizes final output

Helper Function (estimate_noise) [Lines 61-83]:
1. Input Preparation [Lines 62-63]:
   - Ensures noise segment is in column vector format

2. Frame Analysis [Lines 65-71]:
   - Sets up framing parameters
   - Creates Hanning window

3. Spectral Processing [Lines 73-81]:
   - Processes each frame with windowing and FFT
   - Accumulates power spectrum

4. Final Estimation [Line 83]:
   - Averages accumulated spectra for noise profile

#### Main Concepts:
- Short-time Fourier analysis (Lines 32, 76)
- Spectral subtraction for noise reduction (Line 40)
- Overlap-add synthesis (Lines 54-58)
- Window functions for spectral analysis (Lines 24, 29)
- Phase preservation in signal processing (Lines 33, 43)


- **Fourier Transform Plot**
    
```matlab
    function plotfft(audio,fs)
      % Take the FFT of the audio signal
      n = length(audio);  % Number of samples
      f = (0:n-1)*(fs/n);  % Frequency range
      y = fft(audio);  % Compute the FFT

      % Plot the magnitude of the FFT
      plot(f, abs(y));
      title('Magnitude of FFT of Audio Signal');
      xlabel('Frequency (Hz)');
      ylabel('Magnitude');
      xlim([0 fs/2]);  % Plot up to the Nyquist frequency
      grid on;
   end
```

- **Effects of Frequency on Audio**
    
```matlab
    function filteredAudio = filter(audio,fs)
    % Plays the audio signal
    sound(audio, fs);

    % Perform FFT on the audio signal
    Y = fft(audio);

    % Get the length of the audio and calculate frequency vector
    N = length(audio);
    freq = (0:N-1) * (fs / N);

    % Create a mask for the desired frequency range
    mask = (freq <= 7000) & (freq >= 100);
    % For negative frequencies (second half of the FFT)
    mask = mask | (freq >= fs - 7000 & freq <= fs - 100);

    figure;
    % Plot the frequency domain of the original audio
    subplot(2,1,1);
    plot(freq, abs(Y));
    title('Frequency Domain of Original Audio');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 fs/2]);  % Plot up to the Nyquist frequency
    grid on;

    % Zero out frequencies outside the desired range
    Y(~mask) = 0;

    % Plot the frequency domain of Y
    subplot(2,1,2);
    plot(freq, abs(fft(filteredAudio)));
    title('Frequency Domain of Filtered Audio');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 fs/2]);  % Plot up to the Nyquist frequency
    grid on;

    % Inverse FFT to get back to time domain
    filteredAudio = real(ifft(Y));

    % Normalize the audio to prevent clipping
    filteredAudio = filteredAudio / max(abs(filteredAudio));

    % Play the filtered audio
    sound(filteredAudio, fs);
    end
```

The above function helps us to find frequency range for which the audio remains unchanged by listening as well as plotting the frequency domain of the original and filtered audio.

- **Energy Plots (Variation of Energy with Time)**

```matlab
   function energyPlot(audio, fs)
    % Window size (adjust as needed)
    windowSize = round(0.01 * fs);  % 10 ms window
    
    % Overlap between windows
    % overlap = round(windowSize / 2);
    overlap = 0;
    
    % Initialize energy array
    energyOverTime = zeros(1, ceil(length(audio) / (windowSize - overlap)));
    
    % Calculate energy for each window
    for i = 1:length(energyOverTime)
        % Calculate start and end of current window
        startIdx = 1 + (i-1) * (windowSize - overlap);
        endIdx = min(startIdx + windowSize - 1, length(audio));
        
        % Extract window
        window = audio(startIdx:endIdx);
        
        % Calculate instantaneous energy (squared amplitude)
        energyOverTime(i) = sum(window.^2);
    end
    
    % Create time axis
    timeAxis = ((0:length(energyOverTime)-1) * (windowSize - overlap)) / fs;
    
    % Plot energy over time
    figure;
    plot(timeAxis, energyOverTime);
    title('Energy vs Time');
    xlabel('Time (seconds)');
    ylabel('Energy (Squared Amplitude)');
    grid on;
   end
```