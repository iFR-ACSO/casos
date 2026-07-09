rootDir = pwd;
fprintf('Analyzing files from root directory: %s\n\n', rootDir);

% Set up paths for configuration
configFile = fullfile(rootDir, 'casos-coding-guidelines', 'codeAnalyzerConfiguration.json');

% Load and validate configuration
config = struct();
if exist(configFile, 'file')
    fprintf('Loading configuration from: %s\n', configFile);
    try
        % Read config first for custom reporting
        fid = fopen(configFile, 'r');
        configJson = fread(fid, '*char')';
        fclose(fid);
        config = jsondecode(configJson);

        % For MATLAB's code analysis to use the config, we need to either:
        % Option 1: Copy the config to the expected location
        resourcesDir = fullfile(rootDir, 'resources');
        if ~exist(resourcesDir, 'dir')
            mkdir(resourcesDir);
        end
        destConfigFile = fullfile(resourcesDir, 'codeAnalyzerConfiguration.json');
        copyfile(configFile, destConfigFile);

        % Now validate using the resources location
        matlab.codeanalysis.validateConfiguration(rootDir);

        % Refresh to ensure it's used by codeIssues
        matlab.codeanalysis.refreshConfiguration;

        fprintf('Configuration loaded and copied to resources folder successfully\n\n');

        % Display configuration summary
        if isfield(config, 'excludePaths')
            fprintf('Excluding paths: %s\n', strjoin(config.excludePaths, ', '));
        end
        if isfield(config, 'includePaths')
            fprintf('Including paths: %s\n', strjoin(config.includePaths, ', '));
        end
    catch ME
        fprintf('Warning: Configuration issue: %s\n', ME.message);
        fprintf('Using default settings\n\n');
    end
else
    fprintf('No configuration file found at: %s\n', configFile);
    fprintf('Using default settings\n\n');
end

% Run codeIssues - automatically respects the configuration!
fprintf('Running code analysis...\n');
allIssues = codeIssues(rootDir);
fprintf('✓ Analysis complete\n\n');

% Display codeIssues object properties
fprintf('Code Issues Object Properties:\n');
fprintf('  Date: %s\n', allIssues.Date);
fprintf('  Release: %s\n', allIssues.Release);
fprintf('  Configuration: %s\n', allIssues.CodeAnalyzerConfiguration);
fprintf('  Files analyzed: %d\n', length(allIssues.Files));
fprintf('  Active issues: %d\n', height(allIssues.Issues));
fprintf('  Suppressed issues: %d\n\n', height(allIssues.SuppressedIssues));

% Get issues tables
issuesTable = allIssues.Issues;
suppressedTable = allIssues.SuppressedIssues;

% Calculate statistics from the Issues table
if ~isempty(issuesTable)
    % Get unique files with issues
    [uniqueFiles, ~, fileIdx] = unique(issuesTable.Location);
    filesWithIssues = length(uniqueFiles);
    totalIssues = height(issuesTable);

    % Count by severity
    if ~isempty(issuesTable) && any(strcmp(issuesTable.Properties.VariableNames, 'Severity'))
        severityCounts = groupsummary(issuesTable, 'Severity');
    else
        severityCounts = table();
    end

    % Count by check ID (type of issue)
    if any(strcmp(issuesTable.Properties.VariableNames, 'CheckID'))
        checkIDCounts = groupsummary(issuesTable, 'CheckID');
    else
        checkIDCounts = table();
    end
else
    filesWithIssues = 0;
    totalIssues = 0;
    severityCounts = table();
    checkIDCounts = table();
end

% Calculate suppressed issues statistics
if ~isempty(suppressedTable)
    filesWithSuppressions = length(unique(suppressedTable.Location));
    totalSuppressed = height(suppressedTable);
else
    filesWithSuppressions = 0;
    totalSuppressed = 0;
end

% Total MATLAB files
mFiles = dir('**/*.m');
totalFiles = length(mFiles);
filesClean = totalFiles - filesWithIssues;

% Generate Markdown report
reportFile = 'code-analysis-report.md';
fid = fopen(reportFile, 'w');

% Write report header
fprintf(fid, '### MATLAB Code Analysis Report\n\n');
fprintf(fid, '**Generated:** %s  \n', datestr(now));
fprintf(fid, '**MATLAB Release:** %s  \n', allIssues.Release);

% Configuration section
fprintf(fid, '### Analysis Configuration\n\n');
if ~isempty(fieldnames(config))
    fprintf(fid, '```json\n');
    % Pretty print JSON
    configText = jsonencode(config, 'PrettyPrint', true);
    % Replace escaped slashes for readability
    configText = strrep(configText, '\/', '/');
    fprintf(fid, '%s\n', configText);
    fprintf(fid, '```\n\n');
else
    fprintf(fid, '*Default configuration (no config file found)*\n\n');
end

% Summary section
fprintf(fid, '### Summary\n\n');
fprintf(fid, '|        | Value |\n');
fprintf(fid, '|--------|-------|\n');
fprintf(fid, '| **Total Files** | %d |\n', totalFiles);
fprintf(fid, '| **Files Analyzed** | %d |\n', length(allIssues.Files));
fprintf(fid, '| **Clean Files** | %d |\n', filesClean);
fprintf(fid, '| **Files with Issues** | %d |\n', filesWithIssues);
fprintf(fid, '| **Total Active Issues** | %d |\n', totalIssues);
fprintf(fid, '| **Files with Suppressions** | %d |\n', filesWithSuppressions);
fprintf(fid, '| **Total Suppressed Issues** | %d |\n', totalSuppressed);

% Add footer
fprintf(fid, '\n---\n');
fprintf(fid, '*Report generated by MATLAB Code Analyzer*\n');

fclose(fid);

% Generate JSON report with full details
jsonFile = 'code-analysis-report.json';
jsonData = struct();
jsonData.timestamp = datestr(now);
jsonData.matlabRelease = allIssues.Release;
jsonData.configurationStatus = allIssues.CodeAnalyzerConfiguration;
jsonData.stats = struct();
jsonData.stats.totalFiles = totalFiles;
jsonData.stats.filesAnalyzed = length(allIssues.Files);
jsonData.stats.cleanFiles = filesClean;
jsonData.stats.filesWithIssues = filesWithIssues;
jsonData.stats.totalIssues = totalIssues;
jsonData.stats.filesWithSuppressions = filesWithSuppressions;
jsonData.stats.totalSuppressed = totalSuppressed;

% Add configuration
jsonData.configuration = config;

% Add issues (converted from table to struct array)
if totalIssues > 0
    % Convert table to struct array
    issueStruct = table2struct(issuesTable);
    jsonData.issues = issueStruct;
else
    jsonData.issues = [];
end

% Add suppressed issues
if ~isempty(suppressedTable)
    jsonData.suppressedIssues = table2struct(suppressedTable);
else
    jsonData.suppressedIssues = [];
end

% Write JSON file
fid2 = fopen(jsonFile, 'w');
jsonText = jsonencode(jsonData, 'PrettyPrint', true);
jsonText = strrep(jsonText, '\/', '/'); % Clean up escaped slashes
fprintf(fid2, '%s', jsonText);
fclose(fid2);

% Print summary to console
fprintf('\n');
fprintf('========================================\n');
fprintf('  MATLAB CODE ANALYSIS COMPLETE\n');
fprintf('========================================\n');
fprintf(' Date: %s\n', allIssues.Date);
fprintf(' MATLAB Release: %s\n', allIssues.Release);
fprintf(' Configuration: %s\n', allIssues.CodeAnalyzerConfiguration);
fprintf('\n');
fprintf(' Files: %d total, %d analyzed\n', totalFiles, length(allIssues.Files));
fprintf(' Clean files: %d\n', filesClean);
fprintf(' Files with issues: %d\n', filesWithIssues);
fprintf(' Total active issues: %d\n', totalIssues);
fprintf(' Total suppressed issues: %d\n', totalSuppressed);
fprintf('\n');
fprintf(' Reports generated:\n');
fprintf('   %s (Markdown - detailed)\n', reportFile);
fprintf('   %s (JSON - data)\n', jsonFile);
fprintf('========================================\n');