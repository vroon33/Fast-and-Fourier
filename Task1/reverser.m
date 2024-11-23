% Read the audio file
[y, Fs] = audioread('F7.wav');

% Reverse the audio data
y_reversed = flipud(y);

% Save the reversed audio
audiowrite('F7rev.wav', y_reversed, Fs);