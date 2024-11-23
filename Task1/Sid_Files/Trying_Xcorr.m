% PLotting the time domain signal and finding similarity using xcorr

% Load the data
% Checks for the best bird match

clc, clearvars;

% Read the audio file
[y1, Fs1] = audioread('bird1.wav');
[y2, Fs2] = audioread('bird2.wav');
[y3, Fs3] = audioread('bird3.wav');
[y4, Fs4] = audioread('bird3.wav');

v1 = xcorr(y4, y1);
v2 = xcorr(y4, y2);
v3 = xcorr(y4, y3);

score1 = sum(v1);
score2 = sum(v2);
score3 = sum(v3);

% Display the scores
fprintf('Bird 1 Score: %.2f\n', score1);
fprintf('Bird 2 Score: %.2f\n', score2);
fprintf('Bird 3 Score: %.2f\n', score3);