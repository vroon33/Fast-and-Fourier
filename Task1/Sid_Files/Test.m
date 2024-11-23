% Two example vectors
x = [1; 2; 3; 4; 5];
y = [0; 10; 6; 5; -4];

% Compute similarity score
similarityScore = xcorr(x', y')

% Display result
% disp(['Similarity Score: ', num2str(similarityScore)]);
