% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

% Get the root directory
rootDir = pwd;

% Add all necessary paths
addpath(genpath(fullfile(rootDir, 'tests')));
addpath(genpath(rootDir));

disp(rootDir)

% Run tests from tests/core directory using full path
cd(fullfile(rootDir, 'tests', 'solvers'));

% Verify reference file exists
refFile = fullfile(rootDir, 'tests', 'references', 'reference_values.mat');
if exist(refFile, 'file')
    disp('Reference file found!');
else
    error('Reference file not found at: %s', refFile);
end

% Run tests silently (capture results without displaying everything)
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.DiagnosticsRecordingPlugin;

% Create test suite from current directory
suite = matlab.unittest.TestSuite.fromFolder('.');

% Create a silent runner
runner = TestRunner.withNoPlugins();

% Add the DiagnosticsRecordingPlugin (this is the key step!)
runner.addPlugin(DiagnosticsRecordingPlugin());

% Run tests and capture results
results = runner.run(suite);

% ============================================================
% CREATE GITHUB STEP SUMMARY
% ============================================================

%% Get GITHUB_STEP_SUMMARY environment variable
summaryFile = getenv('GITHUB_STEP_SUMMARY');

if isempty(summaryFile)
    warning('GITHUB_STEP_SUMMARY environment variable not found. Running in local mode.');
    summaryFile = 'test_summary_solvers.md';  % Fallback for local testing
end

% Open file for writing (append mode)
fid = fopen(summaryFile, 'a');
if fid == -1
    error('Could not open GITHUB_STEP_SUMMARY file: %s', summaryFile);
end

exitCode = 0;

try
    % Filter results
    failedResults = results([results.Failed]);
    passedResults = results(~[results.Failed]);

    totalTests = numel(results);
    totalFailed = sum([results.Failed]);
    totalPassed = totalTests - totalFailed;
    passRate = (totalPassed / totalTests) * 100;

    if totalFailed > 0
        exitCode = 1;
    end

    % ============================================================
    % WRITE HEADER
    % ============================================================
    fprintf(fid, '#   MATLAB Test Results: Solvers\n\n');
    fprintf(fid, '##  Summary\n\n');
    fprintf(fid, '| Metric | **Total Tests** | **✅ Passed** | **❌ Failed** | **Pass Rate** |\n');
    fprintf(fid, '|--------|-----------------|---------------|---------------|----------------|\n');
    fprintf(fid, '| Value | %d | %d | %d | %.1f%% |\n', totalTests, totalPassed, totalFailed, passRate);

    % ============================================================
    % WRITE FAILED TESTS TABLE 
    % ============================================================

    fprintf(fid, '## Tests Status\n\n');

    % get all of the classes
    flag = zeros(length(results), 1);
    id_plus = 0;
    str_collection = {};
    for i = 1:length(results)
        str = results(i).Name;
        str_split = strsplit(str, '/');
        str_main = str_split{1};
        if ~ismember(str_main, str_collection)
            str_collection{end+1} = str_main;
            id_plus = id_plus+1;
        end
        flag(i)=id_plus;
    end

    % get all of the classes
    flag_failed = [results.Failed];
    flag_passed = [results.Passed];

    % loop over each collection
    for i = 1:length(str_collection)
        idx = (flag == i); % get index of the tests from current test

        total_fails = sum(flag_failed(idx));
        total_succs = sum(flag_passed(idx));

        fprintf(fid,'<details>\n');
        if total_fails > 0
            fprintf(fid,'<summary><code>%s (✅ %d/ ❌ %d )</code></summary>\n', ...
                str_collection{i}, total_succs, total_fails);
        else
            fprintf(fid,'<summary><code>%s (✅ %d)</code></summary>\n', ...
                str_collection{i}, total_succs);
        end
        if total_fails > 0
            fprintf(fid,'<table>\n');
            fprintf(fid,'<thead>\n');
            fprintf(fid,'<tr>\n');
            fprintf(fid,'<th>Test Name</th>\n');
            fprintf(fid,'<th>Failure Reason</th>\n');
            fprintf(fid,'</tr>\n');
            fprintf(fid,'</thead>\n');

            fprintf(fid,'<tbody>\n');
            tempResults = results(idx & flag_failed');

            for j = 1:length(tempResults)
                testName = tempResults(j).Name;

                % Extract failure reason
                failureReason = 'No details available';
                if ~isempty(tempResults(j).Details) && ...
                        isfield(tempResults(j).Details, 'DiagnosticRecord')
                    records = tempResults(j).Details.DiagnosticRecord;
                    if ~isempty(records)
                        % Get the first failure message and clean it
                        try
                            failureReason = records(1).FrameworkDiagnosticResults.DiagnosticText;
                        catch
                            failureReason = records.Exception.message;
                        end
                        % Replace newlines and pipes for Markdown table
                        failureReason = strrep(failureReason, newline, ' ');
                        failureReason = strrep(failureReason, '|', '\\|');
                        % Truncate if too long
                        if length(failureReason) > 100
                            failureReason = [failureReason(1:100) '...'];
                        end
                    end
                end
                % Escape test name for Markdown
                testName = strrep(testName, '|', '\\|');

                fprintf(fid,'<tr>\n');
                fprintf(fid,'<td><code>%s</code></td>\n', testName);
                fprintf(fid,'<td>%s</td>\n', failureReason);
                fprintf(fid,'</tr>\n');
            end
            fprintf(fid,'</tbody>\n');
            fprintf(fid,'</table>\n');
        end
        fprintf(fid,'</details>\n');
    end

    fprintf(fid, '\n');
catch ME
    exitCode = 1;
    fclose(fid);
    rethrow(ME);
end
