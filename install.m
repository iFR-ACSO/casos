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

%% check for casos

% Derive the path for the casos root folder
fullPath = mfilename('fullpath');

% Extract path for the root folder
[folderPath,~,~]=fileparts(fullPath);

% check for all toolboxes on the path 
result = path(); 

% Check if the casos root folder is in the MATLAB path
CasosIsInPath = any(contains(path, folderPath));

% check if if casos is already in matlab path
if CasosIsInPath 

    % fprintf('casos root folder is already on Matlab path.\n')
else
    fprintf('casos root folder added on Matlab path.\n')
    addpath(folderPath);
    savepath(); 
end

%% check for casadi

% check if casadi is installed
casadiStr = which('casadiMEX.mexw64'); 

if isempty(casadiStr)
    error('CasADi is not installed! Please install first')
else
    % fprintf('CasADi is installed\n')

    % check if it is also on the path
    CasadiIsInPath = contains(path, 'casadi');
    if  CasadiIsInPath
        % fprintf('CasADi is also in matlab path.\n')
    else
        warning('CasADi is not in matlab path. It will be added now!')
        % find mosek main folder
        casadiPath = erase(casadiStr,"\casadiMEX.mexw64");

        % add mosek to path and save path
        addpath(genpath(casadiPath))
        savepath();
    end
end

%% check for solver; check for the function that actually calls the optimizer
% mosek
mosekStr = which('mosekopt'); 

if isempty(mosekStr)
    warning('Mosek is not installed')
else
    % fprintf('Mosek is installed\n')

    % check if it is also on the path
    MosekIsInPath = contains(path, 'mosek');
    if  MosekIsInPath
        % fprintf('Mosek is also in matlab path.\n')
    else
        warning('Mosek is not in matlab path. It will be added now!')
        % find mosek main folder
        MosekPath = erase(mosekStr,"\10.1\toolbox\r2017a\mosekopt.mexw64");

        % add mosek to path and save path
        addpath(genpath(MosekPath))
        savepath();
    end
end


% sedumi
sedumiStr = which('sedumi');

if isempty(sedumiStr)
    warning('Sedumi is not installed')
else
    % fprintf('Sedumi is installed\n')

    % check if it is also on the path
    SedumiIsInPath = contains(path, 'sedumi');

    if  SedumiIsInPath
        % fprintf('Mosek is also in matlab path.')
    else
        warning('Mosek is not in matlab path. It will be added now!')
        % find sedumi folder
        SedumiPath = erase(sedumiStr,"sedumi.m");
        % add sedumi to path and save path
        addpath(genpath(SedumiPath))
        savepath();
    end
end


if isempty(sedumiStr) && isempty(mosekStr)
    error('No solver is currently installed! Please install one first')
end

if CasosIsInPath && CasadiIsInPath && (MosekIsInPath || SedumiIsInPath)
    fprintf('CaSos is ready for use!\n')
end