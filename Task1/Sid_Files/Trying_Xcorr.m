clc, clearvars;

[y1, Fs1] = audioread('bird1.wav');
[y2, Fs2] = audioread('bird2.wav');
[y3, Fs3] = audioread('bird3.wav');
[y4, Fs4] = audioread('F7.wav');

v1 = max(xcorr(y4, y1));
v2 = max(xcorr(y4, y2));
v3 = max(xcorr(y4, y3));

fprintf('Bird 1 Score: %.2f\n', v1);
fprintf('Bird 2 Score: %.2f\n', v2);
fprintf('Bird 3 Score: %.2f\n', v3);