function test_runner(silentMode)
    % test_runner- Run test suite with optional silent mode. 
    % Usage: 
    %   test_class();       % Default: verbose mode 
    %   test_class(true);   % Silent mode
    %   test_class(false);  % Verbose mode

    % clean terminal
    clc
    
    % set warning messages to off
    warning('off', 'MATLAB:rankDeficientMatrix');

    % Default to false if no argument is provided
    if nargin < 1
        silentMode = false; 
    end

    % Create a test runner (verbose or silent mode)
    if silentMode
        runner = matlab.unittest.TestRunner.withNoPlugins;
    else
        runner = matlab.unittest.TestRunner.withTextOutput;
    end

    % add to path the folder with tests
    activeFile = matlab.desktop.editor.getActiveFilename;

    % Extract the directory part of the full path
    [activeDir, ~, ~] = fileparts(activeFile);

    % List all entries at the same level as the script
    entries = dir(activeDir);

    % Filter out only directories
    isFolder = [entries.isdir];
    folders = entries(isFolder);

    % Remove '.' and '..' directories
    folders = folders(~ismember({folders.name}, {'.', '..'}));

    % Remove directories that start with '+'
    folders = folders(~startsWith({folders.name}, '+'));

    % Create a list of full paths for each folder
    folderPaths = fullfile(activeDir, {folders.name});

    % Add all folders to the MATLAB path
    addpath(folderPaths{:});
    
    % Display the directory of the file
    disp(['The script is located in: ', activeDir]);
    
    % List all .m files that start with 'test_' in the specified directory
    fileNames = {};
    for i = 1:length(folderPaths)
        filePattern = fullfile(folderPaths{i}, 'test_*.m');
        mFiles = dir(filePattern);
    
        % Initialize a cell array to hold the file names
        fileNames = [fileNames , {mFiles.name}];
    end

    % Display the total number of test files
    disp(['Tests located in: ' folderPaths{i}]);
    disp(['Total number of tests: ', num2str(length(fileNames))]);

    % Create a suite of tests.
    unitTest_suite = testsuite(fileNames);

    % run the tests
    results_suite_basicMath = runner.run(unitTest_suite);

    % Display table with status of all tests
    % disp(results_suite_basicMath.table)

    % Display the results
    statusMessages = {'100% success! Time for a victory dance!', ...
                      'Some tests failed. Keep calm and debug on.'};
    if  all([results_suite_basicMath.Passed])
        fprintf('<strong> STATUS </strong>: %s\n', statusMessages{1});
    else
        fprintf('\n<strong> STATUS </strong>: %s\n', statusMessages{2});
    end
    
end


