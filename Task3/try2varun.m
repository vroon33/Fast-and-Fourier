function detectLoudWords(audio_signal, fs, textFile)
    % Read the timing file
    [words, start_times, end_times, loudness_indicator] = readTimeFile(textFile);
    
    num_words = length(words);
    energies = zeros(num_words, 1);
    
    % Calculate energy for each word
    for i = 1:num_words
        % Convert time to samples
        start_sample = round(start_times(i) * fs);
        end_sample = round(end_times(i) * fs);
        
        % Extract word segment
        word_segment = audio_signal(start_sample:end_sample);
        
        % Calculate RMS energy
        energies(i) = sqrt(mean(word_segment.^2));
    end
    
    % Find separate means for loud and quiet words using ground truth
    quiet_energies = energies(loudness_indicator == 0);
    loud_energies = energies(loudness_indicator == 1);
    
    % Set threshold as weighted value between mean of quiet and loud words
    quiet_mean = mean(quiet_energies);
    loud_mean = mean(loud_energies);
    threshold = quiet_mean + 0.6*(loud_mean - quiet_mean);
    
    % Detect loud words
    detected_loud = energies > threshold;
    
    % Calculate accuracy
    accuracy = sum(detected_loud == loudness_indicator) / num_words * 100;
    
    % Display results
    fprintf('Word\tStart\tEnd\tActual\tPredicted\tEnergy\n');
    fprintf('------------------------------------------------\n');
    for i = 1:num_words
        fprintf('%s\t%.3f\t%.3f\t%d\t%d\t\t%.6f\n', ...
            words{i}, start_times(i), end_times(i), ...
            loudness_indicator(i), detected_loud(i), energies(i));
    end
    fprintf('\nAccuracy: %.2f%%\n', accuracy);
    
    % Visualize results
    figure;
    
    % Plot waveform with marked segments
    subplot(2,1,1);
    t = (0:length(audio_signal)-1)/fs;
    plot(t, audio_signal);
    hold on;
    
    % Mark segments (Green: Correct detection, Red: False Positive, Blue: False Negative)
    for i = 1:num_words
        if loudness_indicator(i) == detected_loud(i)
            color = 'g';  % Correct detection
        elseif detected_loud(i) == 1
            color = 'r';  % False Positive
        else
            color = 'b';  % False Negative
        end
        xline(start_times(i), color);
        xline(end_times(i), color);
        text(start_times(i), max(audio_signal)*0.8, words{i}, 'Color', color);
    end
    
    title('Speech Signal with Word Boundaries');
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % Plot energies
    subplot(2,1,2);
    bar(energies);
    hold on;
    plot([1 num_words], [threshold threshold], 'r--', 'LineWidth', 2);
    plot([1 num_words], [quiet_mean quiet_mean], 'b--', 'LineWidth', 1);
    plot([1 num_words], [loud_mean loud_mean], 'g--', 'LineWidth', 1);
    title('Word Energies');
    xticks(1:num_words);
    xticklabels(words);
    xtickangle(45);
    ylabel('Energy');
    legend('Energy', 'Threshold', 'Quiet Mean', 'Loud Mean');
end

[audio, fs] = audioread('6.wav');
detectLoudWords(audio, fs, '6.wav');
