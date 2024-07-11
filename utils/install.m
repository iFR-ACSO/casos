% -------------------------------------------------------------------------
%
% Short Description: Script that install casos by adding it to the matlab
%                    path and checks if casadi is installed. If either 
%                    casos or casadi are not on the path they are added to 
%                    it. Also it is checked which solvers are available and
%                    if they are on the current path.
%
% Date: 04/02/2024
%
% -------------------------------------------------------------------------

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
casadiStr = which('casadiMEX.mexw64');

casadiVersionCorrect = true;

if isempty(casadiStr) | ~exist('casadi.SX', 'class')
    warning(' -- CasADi is not installed! Please install first --\n')
    casadiInstalled = false;
else
    casadiInstalled = true;
    fprintf(' -- CasADi is installed --\n')
    fprintf(' -- Checking CasADi version --\n')
    % check if casadi.DM.sparsity_cast is available
    casadi.DM; % load
    casadiVersionStr = which('sparsity_cast');
    if isempty(casadiVersionStr)
        warning('CasADi version is not up to date')
        casadiVersionCorrect = false;
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
installationOverview = table(software, installed)

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

