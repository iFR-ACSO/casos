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
scriptPath = mfilename('fullpath');

% grandparent folder
casosPath = fileparts(fileparts(scriptPath));

% check if casos.PS class exists
casosExists = exist('casos.PS', 'class');

if ~casosExists
    warning( ...
        "\n --- casos.PS class not available. Adding casos root folder to path ---\n")
    
    % check for casos root folder
    if ~endsWith(casosPath, 'casos')
        error("Could not find casos root folder. Please add to path manually\n")
    end

    % add casos to path
    fprintf("--- Adding casos root folder to path ---\n")
    addpath(casosPath);
    savepath();

    % check if adding Casos to path worked
    if ~exist('casos.PS', 'class')
        error("Unable to add Casos to path. Please add manually\n");
    end

    fprintf("--- Casos was added to path ---\n")
else
    fprintf("\n--- Casos already installed ---\n")
end


%% check for casadi
% check if casadi is installed
casadiStr = which('casadiMEX.mexw64'); 

if isempty(casadiStr) | ~exist('casadi.SX', 'class')
    warning('--- CasADi is not installed! Please install first ---\n')
    casadiInstalled = false;
else
    fprintf('--- CasADi is installed ---\n')
    casadiInstalled = true;
end

% 
% %% check for solver; check for the function that actually calls the optimizer
% % mosek
% mosekStr = which('mosekopt'); 
% 
% if isempty(mosekStr)
%     warning('Mosek is not installed')
% else
%     % fprintf('Mosek is installed\n')
% 
%     % check if it is also on the path
%     MosekIsInPath = contains(path, 'mosek');
%     if  MosekIsInPath
%         % fprintf('Mosek is also in matlab path.\n')
%     else
%         warning('Mosek is not in matlab path. It will be added now!')
%         % find mosek main folder
%         MosekPath = erase(mosekStr,"\10.1\toolbox\r2017a\mosekopt.mexw64");
% 
%         % add mosek to path and save path
%         addpath(genpath(MosekPath))
%         savepath();
%     end
% end
% 
% 
% % sedumi
% sedumiStr = which('sedumi');
% 
% if isempty(sedumiStr)
%     warning('Sedumi is not installed')
% else
%     % fprintf('Sedumi is installed\n')
% 
%     % check if it is also on the path
%     SedumiIsInPath = contains(path, 'sedumi');
% 
%     if  SedumiIsInPath
%         % fprintf('Mosek is also in matlab path.')
%     else
%         warning('Mosek is not in matlab path. It will be added now!')
%         % find sedumi folder
%         SedumiPath = erase(sedumiStr,"sedumi.m");
%         % add sedumi to path and save path
%         addpath(genpath(SedumiPath))
%         savepath();
%     end
% end
% 
% 
% if isempty(sedumiStr) && isempty(mosekStr)
%     error('No solver is currently installed! Please install one first')
% end
% 
% if CasosIsInPath && CasadiIsInPath && (MosekIsInPath || SedumiIsInPath)
%     fprintf('CaSos is ready for use!\n')
% end

% print table
software = ["Casos"; "Casadi"; "Mosek"; "Sedumi"];
installed = [true; true; true; true];

installationOverview = table(software, installed)