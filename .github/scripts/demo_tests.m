% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

% Add all necessary paths to MATLAB search path
addpath(genpath(pwd));

% Navigate to demo folder containing all example scripts
cd('demo')

% Find all MATLAB files in the demo directory
mfiles = dir('*.m');

% Initialize arrays to track results
success_files = {};
failed_files = {};
error_messages = {};
flag_failed = zeros(length(mfiles),1);

% Save initial workspace variables that we want to preserve
track_vars = {'mfiles', 'success_files', 'failed_files', 'error_messages', 'track_vars', 'i', 'flag_failed'};

runDemo = @(filename) run(filename);

% Execute each demo file individually
% Using try-catch to continue even if a demo fails
for ifiles=1:length(mfiles)
    fprintf('\n# ========================================= #\n')
    fprintf('  %s\n', mfiles(ifiles).name)
    fprintf('# ========================================= #\n')

    try
        % Run the demo script
        run(mfiles(ifiles).name);
        fprintf('%s: successfull!', mfiles(ifiles).name)
        success_files{end+1} = mfiles(ifiles).name;
    catch ME
        % Capture error details but continue to next demo
        disp(ME);
        warning('File %s encountered an error!', mfiles(ifiles).name)
        flag_failed(ifiles) = 1;
        failed_files{end+1} = mfiles(ifiles).name;
        error_messages{end+1} = sprintf('%s: %s', mfiles(ifiles).name, ME.message);
    end
    % Clear workspace but preserve tracking variables
    % Get list of all variables in current workspace
    all_vars = who;
    % Remove tracking variables from the list
    vars_to_clear = setdiff(all_vars, track_vars);
    % Clear only the non-tracking variables
    clear(vars_to_clear{:});
end

%% ============================================
% GENERATE GITHUB STEP SUMMARY
% ============================================
summaryFile = getenv('GITHUB_STEP_SUMMARY');

if isempty(summaryFile)
    warning('GITHUB_STEP_SUMMARY environment variable not found. Running in local mode.');
    summaryFile = 'test_summary_demos.md';  % Fallback for local testing
end

% Open file for writing (append mode)
fid = fopen(summaryFile, 'a');
if fid == -1
    error('Could not open GITHUB_STEP_SUMMARY file: %s', summaryFile);
end

exitCode = 0;

try
    fprintf(fid, '#   MATLAB Test Results: Demos\n\n');
    fprintf(fid, '##  Summary\n\n');
    fprintf(fid, '| Metric | **Total Tests** | **✅ Passed** | **❌ Failed** | **Pass Rate** |\n');
    fprintf(fid, '|--------|-----------------|---------------|---------------|----------------|\n');
    fprintf(fid, '| Value | %d | %d | %d | %.1f%% |\n', length(mfiles), length(success_files), ...
        length(failed_files), 100*length(success_files)/length(mfiles));

    if ~isempty(failed_files)
        exitCode = 1;
    end

    % ============================================================
    % WRITE FAILED TESTS TABLE (ONLY IF THERE ARE FAILURES)
    % ============================================================
    fprintf(fid, '## Tests Status\n\n');

    idx_error_message = 0;
    for i=1:length(mfiles)

        fprintf(fid,'<details>\n');

        if flag_failed(i)==1
            fprintf(fid,'<summary><code>%s (❌)</code></summary>\n', ...
                mfiles(i).name);
        else
            fprintf(fid,'<summary><code>%s (✅)</code></summary>\n', ...
                mfiles(i).name);
        end


        if flag_failed(i)==1
            fprintf(fid,'<table>\n');
            fprintf(fid,'<thead>\n');
            fprintf(fid,'<tr>\n');
            fprintf(fid,'<th>Failure Reason</th>\n');
            fprintf(fid,'</tr>\n');
            fprintf(fid,'</thead>\n');
        
            fprintf(fid,'<tbody>\n');
        
            fprintf(fid,'<tr>\n');
    
            idx_error_message = idx_error_message+1;
            failureReason = strrep(error_messages{idx_error_message}, newline, ' ');
            failureReason = strrep(failureReason, '|', '\\|');
            % Truncate if too long
            if length(failureReason) > 500
                failureReason = [failureReason(1:500) '...'];
            end
            fprintf(fid,'<td>%s</td>\n', failureReason);
                    fprintf(fid,'</tr>\n');
            
            
        end
        fprintf(fid,'</tbody>\n');
        fprintf(fid,'</table>\n');
        fprintf(fid,'</details>\n');
        fprintf(fid, '\n');
    end
catch ME
    exitCode = 1;
    fclose(fid);
    rethrow(ME);
end