clc
clear


%% setup a testsuite for basic math operations

% add to path the folder with tests
activeFile = matlab.desktop.editor.getActiveFilename;

% Extract the directory part of the full path
[activeDir, ~, ~] = fileparts(activeFile);
test_folder_dir = [activeDir, '/PolyOperations'];
addpath(test_folder_dir);

% Display the directory of the file
disp(['The script is located in: ', activeDir]);

% List all .m files that start with 'test_' in the specified directory
filePattern = fullfile(test_folder_dir, 'test_*.m');
mFiles = dir(filePattern);

% Create a suite of tests.
suite_basicMath = testsuite({mFiles.name});

% run the tests
results_suite_basicMath = run(suite_basicMath);

if  all([results_suite_basicMath.Passed])
    disp('Passed all test for basic polynomial math operations.')
else
  
end




