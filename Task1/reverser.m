% Read the audio file
[y, Fs] = audioread('bird3.wav');

% Reverse the audio data
y_reversed = flipud(y);

% Save the reversed audio
audiowrite('F7_rev.wav', y_reversed, Fs);