function [words, startTimes, endTimes, isLouder] = readTimeFile(filename)
    % Read the entire file as a string
    fileContent = fileread(filename);
    
    % Split the content by newlines or spaces to get individual elements
    % Using regexp to split by any number of spaces, tabs, or newlines
    elements = regexp(fileContent, '\s+', 'split');
    
    % Remove any empty elements
    elements = elements(~cellfun('isempty', elements));
    
    % Initialize arrays
    numWords = length(elements)/4;  % Each word has 4 components
    words = cell(numWords, 1);
    startTimes = zeros(numWords, 1);
    endTimes = zeros(numWords, 1);
    isLouder = zeros(numWords, 1);
    
    % Parse the elements
    for i = 1:numWords
        baseIdx = (i-1)*4 + 1;
        words{i} = elements{baseIdx};
        startTimes(i) = str2double(elements{baseIdx + 1});
        endTimes(i) = str2double(elements{baseIdx + 2});
        isLouder(i) = str2double(elements{baseIdx + 3});
    end
    
    % Display parsed data to verify
    for i = 1:length(words)
        fprintf('%s\t\t%f\t%f\t%d\n', words{i}, startTimes(i), endTimes(i), isLouder(i));
    end
    fprintf('\n');
end