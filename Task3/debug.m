% Print current directory and list contents
currentDir = pwd;
fprintf('Current directory: %s\n', currentDir);

% List contents of current directory
fprintf('\nContents of current directory:\n');
dir

% List contents of parent directory
fprintf('\nContents of parent directory:\n');
parentContents = dir('..');
for i = 1:length(parentContents)
    fprintf('%s\n', parentContents(i).name);
end

% Try to find the audio file using different possible paths
possiblePaths = {
    fullfile('..', 'References and Tasks', 'Project_LouderWordsDetection', 'audios', '2.wav'),
    fullfile('..', 'References_and_Tasks', 'Project_LouderWordsDetection', 'audios', '2.wav'),
    fullfile('..', 'References & Tasks', 'Project_LouderWordsDetection', 'audios', '2.wav'),
    fullfile('..', 'References_Tasks', 'Project_LouderWordsDetection', 'audios', '2.wav')
};

% Try each possible path
fprintf('\nTrying possible paths:\n');
for i = 1:length(possiblePaths)
    fprintf('Trying path: %s\n', possiblePaths{i});
    if exist(possiblePaths{i}, 'file') == 2
        fprintf('Found file at: %s\n', possiblePaths{i});
        audioPath = possiblePaths{i};
        break;
    end
end

% If we found the file, proceed with processing
if exist('audioPath', 'var')
    [audio, fs] = audioread(audioPath);
    fprintf('Successfully read audio file!\n');
    
    % Process the audio as before
    % ... rest of your processing code ...
else
    fprintf('Could not find audio file in any of the expected locations.\n');
    fprintf('Please verify the correct path structure and update the code accordingly.\n');
end