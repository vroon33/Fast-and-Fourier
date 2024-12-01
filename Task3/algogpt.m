% Load the audio file
[audioSignal, fs] = audioread('2.wav'); % Replace with your file
audioSignal = audioSignal(:, 1); % Use mono signal if stereo

% Parameters
frameDuration = 0.02; % 20 ms frame
frameLength = round(frameDuration * fs); % Number of samples per frame
overlap = round(0.5 * frameLength); % 50% overlap
thresholdEnergy = 0.01; % Adjust as per your audio signal
thresholdZCR = 0.15; % Adjust as per your audio signal

% Frame the signal
window = hamming(frameLength); % Hamming window
numFrames = floor((length(audioSignal) - frameLength) / (frameLength - overlap)) + 1;
frames = buffer(audioSignal, frameLength, overlap, 'nodelay') .* window;

% Compute Short-Time Energy (STE)
STE = sum(frames .^ 2, 1);

% Compute Zero-Crossing Rate (ZCR)
ZCR = sum(abs(diff(sign(frames))) > 0, 1) / frameLength;

% Normalize features for thresholding
STE = STE / max(STE);
ZCR = ZCR / max(ZCR);

% Detect speech frames
speechFrames = (STE > thresholdEnergy) & (ZCR < thresholdZCR);

% Convert frames to timestamps
speechStartEnd = [];
for i = 1:length(speechFrames)
    if i == 1 || speechFrames(i) ~= speechFrames(i-1)
        speechStartEnd = [speechStartEnd; i];
    end
end
if mod(length(speechStartEnd), 2) ~= 0
    speechStartEnd = [speechStartEnd; length(speechFrames)];
end

% Convert frame indices to time
speechTimestamps = (speechStartEnd - 1) * (frameLength - overlap) / fs;

% Display timestamps
disp('Speech segments detected at (seconds):');
disp(speechTimestamps);

% Plot signal with detected segments
time = (0:length(audioSignal)-1) / fs;
figure;
plot(time, audioSignal);
hold on;
for i = 1:2:length(speechTimestamps)
    % Ensure indices are within valid range
    startIdx = max(1, round(speechTimestamps(i) * fs));
    endIdx = min(length(audioSignal), round(speechTimestamps(i+1) * fs));
    
    % Extract the segment
    x = (startIdx:endIdx) / fs;
    y = audioSignal(startIdx:endIdx);
    plot(x, y, 'r', 'LineWidth', 2);
end
xlabel('Time (s)');
ylabel('Amplitude');
title('Detected Speech Segments');
legend('Audio Signal', 'Speech Segments');
grid on;
