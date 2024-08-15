% -------------------------------------------------------------------------
% Casos Installation Check
% -------------------------------------------------------------------------
% 
% Purpose:
%   Verifies the installation and configuration of essential toolboxes 
%   (Casos, CasADi) and solvers (Mosek, SeDuMi) required for Casos.
% 
% Usage:
%   Run this script to ensure all dependencies for Casos are correctly 
%   installed before use.


fprintf(['\n\n----------------------------------------------------\n' ...
    '         Checking installed toolboxes\n' ...
    '----------------------------------------------------\n\n'])

%% Check for Casos in path
% get path of this script
fprintf("--- Checking Casos Installation ---\n")
scriptPath = mfilename('fullpath');

% grandparent folder
casosPath = fileparts(fileparts(scriptPath));

% check if casos.PS class exists
casosExists = exist('casos.PS', 'class');

if ~casosExists
    warning( ...
        " -- casos.PS class not available. Adding casos root folder " + ...
        "to path --\n")
    
    % check for casos root folder
    if ~endsWith(casosPath, 'casos')
        error("Could not find casos root folder. Please add to " + ...
            "path manually\n")
    end

    % add casos to path
    fprintf(" -- Adding casos root folder to path --\n")
    addpath(casosPath);
    savepath();

    % check if adding Casos to path worked
    if ~exist('casos.PS', 'class')
        error("Unable to add Casos to path. Please add manually\n");
    end

    fprintf(" -- Casos was added to path ---\n")
else
    fprintf(" -- Casos already installed --\n")
end


%% check for casadi
% check if casadi is installed
fprintf("\n--- Checking CasADi installation ---\n")
casadiStr = which('casadiMEX');

casadiVersionCorrect = true;

if isempty(casadiStr) | ~exist('casadi.SX', 'class')
    warning(' -- CasADi is not installed! Please install first --\n')
    casadiInstalled = false;
else
    casadiInstalled = true;
    fprintf(' -- CasADi is installed --\n')
    
    % check if casadi.DM.sparsity_cast is available
    casadi.DM; % load
    casadiVersionStr = which('sparsity_cast');
    if isempty(casadiVersionStr)
        warning('CasADi version is not up to date')
        casadiVersionCorrect = false;
    else
        fprintf(' -- CasADi version is up to date --\n')
    end
end


%% Check for solvers
% mosek
fprintf("\n--- Checking Mosek Installation ---\n")

mosekStr = which('mosekopt');

mosekInstalled = true;
if isempty(mosekStr)
    warning('Mosek is not installed')
    mosekInstalled = false;
else
    fprintf(' -- Mosek is installed --\n')
end




% SeDuMi
fprintf("\n--- Checking SeDuMi Installation ---\n")

sedumiStr = which('sedumi');

sedumiInstalled = true;
if isempty(sedumiStr)
    warning('SeDuMi is not installed')
    sedumiInstalled = false;
else
    fprintf(' -- SeDuMi is installed --\n')
end

%% Final results
% Print overview table
software = ["Casos"; "Casadi"; "Mosek"; "SeDuMi"];
installed = [true; casadiInstalled; mosekInstalled; sedumiInstalled];
fprintf('\n--- Summary ---\n\n')
disp(table(software, installed))

% error if CasADi is not installed
if ~casadiInstalled
    error('CasADi needs to be installed.')
end

% error if no solver is installed
if ~mosekInstalled && ~sedumiInstalled
    error('No solver is currently installed! Install at least one solver')
end

% error if CasADi version is too old
if ~casadiVersionCorrect
    error(['Your installation of CasADi is too old. Casos might not ' ...
        'run correctly'])
end
fprintf('--- Casos is ready to be used ---\n')

