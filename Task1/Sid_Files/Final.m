% Checks for the best bird match

clc, clearvars;

% Read the audio file
[y, Fs] = audioread('bird1.wav');
[Top_Peak_Frequencies1, Top_Peak_Magnitudes1] = Peak_Processing(y, Fs);

[y, Fs] = audioread('bird2.wav');
[Top_Peak_Frequencies2, Top_Peak_Magnitudes2] = Peak_Processing(y, Fs);

[y, Fs] = audioread('bird3.wav');
[Top_Peak_Frequencies3, Top_Peak_Magnitudes3] = Peak_Processing(y, Fs);

[y, Fs] = audioread('F1.wav'); % Test file
[Top_Peak_Frequencies_Test, Top_Peak_Magnitudes_Test] = Peak_Processing(y, Fs);

% Find the frequency in Reference files closest to the frequency in test file and substract
[diff1, idx1] = findClosestMultiple(Top_Peak_Frequencies_Test, Top_Peak_Frequencies1);
[diff2, idx2] = findClosestMultiple(Top_Peak_Frequencies_Test, Top_Peak_Frequencies2);
[diff3, idx3] = findClosestMultiple(Top_Peak_Frequencies_Test, Top_Peak_Frequencies3);

idx1 = 11 - idx1;
idx2 = 11 - idx2;
idx3 = 11 - idx3;

% Change value to 0 if the value is greater than 100
diff1(diff1 > 100) = 0;
diff2(diff2 > 100) = 0;
diff3(diff3 > 100) = 0;

% Subtract the value from 100 and add to get individual score, Normalized
score1 = sum((100 - diff1) .* idx1) / 55;
score2 = sum((100 - diff2) .* idx2) / 55;
score3 = sum((100 - diff3) .* idx3) / 55;

% Display the scores
fprintf('Bird 1 Score: %.2f%%\n', score1);
fprintf('Bird 2 Score: %.2f%%\n', score2);
fprintf('Bird 3 Score: %.2f%%\n', score3);

function [closest_values, Index] = findClosestMultiple(x, a)
    % Initialize output arrays
    closest_values = zeros(size(x));
    Index = zeros(1, length(x));
    
    % For each element in x
    for i = 1:length(x)
        % Find the absolute differences between current x and all elements in a
        differences = abs(a - x(i));
        
        % Find the minimum difference and its index
        [close, idx] = min(differences);
        
        % Store results
        closest_values(i) = close;
        Index(i) = idx;
    end
end