% Read the audio file
    [y, fs] = audioread('2.wav');
    
    % Convert stereo to mono if necessary
    if size(y, 2) > 1
        y = mean(y, 2);
    end
    
    % Noise filtering
    y = noisefilter(y, fs); 

    % Parameters for analysis
    windowSize = round(0.001 * fs);  % 30ms window
    overlap = round(0.005 * fs);    % 15ms overlap
    threshold = 3;                % Adjust this value based on your needs
    
    % Calculate short-time energy
    energy = zeros(floor((length(y)-windowSize)/overlap), 1);
    for i = 1:length(energy)
        startIdx = (i-1)*overlap + 1;
        endIdx = startIdx + windowSize - 1;
        frame = y(startIdx:endIdx);
        energy(i) = sum(frame.^2);
    end
    
    % Normalize energy
    energy = energy / max(energy);
    
    % Find segments with high energy (emphasized parts)
    emphasized = energy > threshold * mean(energy);
    
    % Convert frame indices to time
    timeAxis = ((0:length(energy)-1) * overlap) / fs;
    
    % Plot results
    figure;
    
    % Plot waveform
    subplot(3,1,1);
    plot((1:length(y))/fs, y);
    title('Audio Waveform');
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % Plot energy
    subplot(3,1,2);
    plot(timeAxis, energy);
    title('Short-time Energy');
    xlabel('Time (s)');
    ylabel('Normalized Energy');
    
    % Plot emphasized segments
    subplot(3,1,3);
    plot(timeAxis, emphasized);
    title('Emphasized Segments (1 = emphasized)');
    xlabel('Time (s)');
    ylabel('Emphasis Detection');
    ylim([-0.1 1.1]);
    
    % Print time stamps of emphasized segments
    emphasisChanges = diff([0; emphasized; 0]);
    startPoints = find(emphasisChanges == 1);
    endPoints = find(emphasisChanges == -1) - 1;
    
    fprintf('\nEmphasized segments found at:\n');
    for i = 1:length(startPoints)
        startTime = timeAxis(startPoints(i));
        endTime = timeAxis(endPoints(i));
        fprintf('%.2f seconds to %.2f seconds\n', startTime, endTime);
    end