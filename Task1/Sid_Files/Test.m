% Calls the Peak_processing function

% Read the audio file
[y, Fs] = audioread('bird1.wav');

% Call the Peak_Processing function
[Top_Peak_Frequencies, Top_Peak_Magnitudes] = Peak_Processing(y, Fs);

% Display the top 10 peaks and their frequencies
fprintf('Top 10 Peaks:\n');
fprintf('Frequency (Hz): %.2f\tMagnitude: %.2f\n', [Top_Peak_Frequencies; Top_Peak_Magnitudes]);