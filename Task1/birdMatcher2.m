function matchBirdSound()
    % Read reference files
    [bird1, fs1] = audioread('bird1.wav');
    [bird2, fs2] = audioread('bird2.wav');
    [bird3, fs3] = audioread('bird3.wav');
    [taskFile, fsTask] = audioread('F5.wav');
    
    % Ensure all files have same sampling rate
    if fs1 ~= fs2 || fs2 ~= fs3 || fs3 ~= fsTask
        error('All audio files must have the same sampling rate');
    end
    
    % Using Correlation
    corr1 = max(xcorr(bird1, taskFile));
    corr2 = max(xcorr(bird2, taskFile));
    corr3 = max(xcorr(bird3, taskFile));
    
    % Normalize correlation values
    corrValues = [corr1, corr2, corr3];
    corrValues = corrValues / max(abs(corrValues));
    
    % Print results
    fprintf('\nCorrelation Method Results:\n');
    fprintf('Bird 1 match: %.2f%%\n', corrValues(1)*100);
    fprintf('Bird 2 match: %.2f%%\n', corrValues(2)*100);
    fprintf('Bird 3 match: %.2f%%\n', corrValues(3)*100);
    
    % Determine best match
    [bestScore, bestIndex] = max(corrValues);
    fprintf('\nBest Match: Bird %d (Score: %.2f%%)\n', ...
        bestIndex, bestScore*100);
    
    % Plot spectrograms for visual comparison
    figure;
    
    % Define spectrogram parameters
    window = hamming(256);
    noverlap = 128;
    nfft = 512;
    
    subplot(2,2,1);
    spectrogram(taskFile, window, noverlap, nfft, fsTask, 'yaxis');
    title('Task File Spectrogram');
    
    subplot(2,2,2);
    spectrogram(bird1, window, noverlap, nfft, fs1, 'yaxis');
    title('Bird 1 Spectrogram');
    
    subplot(2,2,3);
    spectrogram(bird2, window, noverlap, nfft, fs2, 'yaxis');
    title('Bird 2 Spectrogram');
    
    subplot(2,2,4);
    spectrogram(bird3, window, noverlap, nfft, fs3, 'yaxis');
    title('Bird 3 Spectrogram');
    
    % Plot correlation results in a new figure
    figure;
    bar(corrValues*100);
    title('Bird Sound Matching Results');
    xlabel('Bird Reference Number');
    ylabel('Match Percentage (%)');
    ylim([0 100]);
    grid on;
end