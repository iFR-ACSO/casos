% ========================================================================
%
% Name: test_runner
%
% Short Description: Script executing and checking tests.
%                
% Date: 08/18/2023 
%
%
% ========================================================================

% Get all function files in test f
allFiles = dir(fullfile('./+Tests/', '*.m'));

% Initialize a cell array to store the names of functions
functionNames = char();

% Loop through each file and identify functions
for i = 1:numel(allFiles)
    [~, functionName, ext] = fileparts(allFiles(i).name);
    if strcmp(ext, '.m') && ~isempty(functionName) && ~strcmp(functionName, 'test_script')
        % Exclude empty names and the script's own name if it's named 'test_script.m'
        functionNames{end+1} = functionName;
    end
end


% Initialize a cell array to store the names of functions that failed the test
failedFunctions = {};

% Loop through each function file and execute the function
for i = 1:numel(functionNames)
    functionName = ['Tests.' functionNames{i}]; % Remove '.m' extension
    try
        % Execute the function and get the boolean result
        result = feval(functionName);
        
        % Check if the result is false
        if ~result
            failedFunctions{end+1} = functionName;
        end
    catch
        % Catch any errors that occur during function execution
        failedFunctions{end+1} = functionName;
    end
end

% Display the list of failed functions
if isempty(failedFunctions)
    disp('All functions passed the test.');
else
    disp('Failed functions:');
    disp(failedFunctions);
end