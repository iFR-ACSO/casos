function test_class(silentMode)
    % test_class - Run test suite with optional silent mode.
    % Usage:
    %   test_class();       % Default: verbose mode
    %   test_class(true);   % Silent mode
    %   test_class(false);  % Verbose mode

    % clean workspace and terminal
    clc
    clear
    
    % Default to false if no argument is provided
    if nargin < 1
        silentMode = false; 
    end

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
    
    % Initialize a cell array to hold the file names
    fileNames = {mFiles.name};
    
    % Display the total number of test files
    disp(['Total number of tests: ', num2str(length(fileNames))]);
    
    % Create a suite of tests.
    suite_basicMath = testsuite(fileNames);
    
    % Create a test runner with no output plugins (silent mode).
    runner = matlab.unittest.TestRunner.withNoPlugins;
    
    % run the tests
    if silentMode
        results_suite_basicMath = runner.run(suite_basicMath);
    else
        results_suite_basicMath = run(suite_basicMath);
    end
    
    % Display the results
    if  all([results_suite_basicMath.Passed])
        disp('Passed all test for basic polynomial math operations.')
    else
        disp('Failed some test for basic polynomial math operations.')
    end

end
