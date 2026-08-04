% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function package_toolbox(version)
% Create a zip file.
%
% IMPORTANT: The script requires the current folder to be the casos root
% folder!

assert(isfolder("+casos"), "Current folder must be casos root folder.")

% ensure string
version = string(version);

%% Create temporary folder structure
% create the following folder structure:
%
%   casos/ temp/
%   |-- +casos/
%   |       :
%   |-- doc/
%   |       GettingStarted.mlx
%   |-- demo/
%   |       :

if exist("temp","dir")
    % clear temporary folder and structure
    rmdir("temp","s")
end

mkdir("temp")
mkdir("temp/+casos")
% mkdir("temp/doc")
mkdir("temp/demo")

% copy source folder
copyMfiles("+casos","temp/+casos")
% copy examples from demo folder
copyMfiles("demo","temp/demo")
% copy Getting Started guide
% copyfile("build/GettingStarted.mlx","temp/doc/")
% copy license
copyfile("LICENSES/GPL-3.0-only.txt","temp")

%% Generate meta data
% populate CasosMeta class folder
generate_metadata("temp",version);

%% Create zip file
% use temporary folder as root
zip(join(["build/casos" version "matlab.zip"],"-"), [
    "+casos"
    "demo"
    % "doc"
    "GPL-3.0-only.txt"
], "temp")

end
