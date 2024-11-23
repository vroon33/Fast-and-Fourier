function matchBirdSound()
    % Read reference files
    [bird1, fs1] = audioread('bird1.wav');
    [bird2, fs2] = audioread('bird2.wav');
    [bird3, fs3] = audioread('bird3.wav');
    [taskFile, fsTask] = audioread('F8.wav');
    
    % Ensure all files have same sampling rate
    if fs1 ~= fs2 || fs2 ~= fs3 || fs3 ~= fsTask
        error('All audio files must have the same sampling rate');
    end
    
    % Step 1: Time Domain Correlation (Primary Check)
    corr1 = max(xcorr(bird1, taskFile));
    corr2 = max(xcorr(bird2, taskFile));
    corr3 = max(xcorr(bird3, taskFile));
    
    corrValues = [corr1, corr2, corr3];
    corrValues = corrValues / max(abs(corrValues));
    
    % First decision threshold
    if max(corrValues) > 0.7  % Strong correlation found
        [~, bestIndex] = max(corrValues);
        finalScores = corrValues;
    else
        % Step 2: Energy Pattern Analysis
        window = 1000;
        energy_pattern = @(x) sum(reshape(x(1:floor(length(x)/window)*window), window, []).^2);
        
        energy1 = energy_pattern(bird1);
        energy2 = energy_pattern(bird2);
        energy3 = energy_pattern(bird3);
        energyTask = energy_pattern(taskFile);
        
        energyCorr1 = corr(energy1', energyTask');
        energyCorr2 = corr(energy2', energyTask');
        energyCorr3 = corr(energy3', energyTask');
        
        energyValues = [energyCorr1, energyCorr2, energyCorr3];
        energyValues = (energyValues - min(energyValues)) / (max(energyValues) - min(energyValues));
        
        % Combine with original correlation with more weight to correlation
        finalScores = (0.7 * corrValues + 0.3 * energyValues);
    end
    
    % Apply confidence threshold
    confidenceThreshold = 0.3;
    finalScores(finalScores < confidenceThreshold) = 0;
    
    % Normalize final scores
    if max(finalScores) > 0
        finalScores = finalScores / max(finalScores);
    end
    
    % Print results
    fprintf('\nMatching Results:\n');
    fprintf('Bird 1 match: %.2f%%\n', finalScores(1)*100);
    fprintf('Bird 2 match: %.2f%%\n', finalScores(2)*100);
    fprintf('Bird 3 match: %.2f%%\n', finalScores(3)*100);
    
    % Determine best match
    [bestScore, bestIndex] = max(finalScores);
    
    if bestScore > 0
        fprintf('\nBest Match: Bird %d (Score: %.2f%%)\n', ...
            bestIndex, bestScore*100);
    else
        fprintf('\nNo confident match found\n');
    end
    
    % Plot spectrograms
    figure;
    window_spec = hamming(256);
    noverlap = 128;
    nfft_spec = 512;
    
    subplot(2,2,1);
    spectrogram(taskFile, window_spec, noverlap, nfft_spec, fsTask, 'yaxis');
    title('Task File Spectrogram');
    
    subplot(2,2,2);
    spectrogram(bird1, window_spec, noverlap, nfft_spec, fs1, 'yaxis');
    title('Bird 1 Spectrogram');
    
    subplot(2,2,3);
    spectrogram(bird2, window_spec, noverlap, nfft_spec, fs2, 'yaxis');
    title('Bird 2 Spectrogram');
    
    subplot(2,2,4);
    spectrogram(bird3, window_spec, noverlap, nfft_spec, fs3, 'yaxis');
    title('Bird 3 Spectrogram');
    
    % Plot final scores
    figure;
    bar(finalScores*100);
    title('Final Matching Scores');
    xlabel('Bird Reference Number');
    ylabel('Match Percentage (%)');
    ylim([0 100]);
    grid on;
end