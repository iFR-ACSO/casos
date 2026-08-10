% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function generate_metadata(path_to_source, version)
% Create CaΣoS meta data.

% REUSE-IgnoreStart
license_header = "% SPDX-License-Identifier: CC0-1.0";
% REUSE-IgnoreEnd
path_to_meta = join([path_to_source "+casos" "@CasosMeta"], filesep);

%% Create method to query version
% create a method file
lines_version = compose("function v = version(), v = '%s'; end", version);
writelines([license_header lines_version], join([path_to_meta "version.m"], filesep));

%% Create methods to query git commit
% only if inside a git repository
try
    % git repository
    repo = gitrepo;
    last_commit = repo.LastCommit;

    % create method to query git revision
    lines_revision = compose("function v = git_revision(), v = '%s'; end", last_commit.ID);
    writelines([license_header lines_revision], join([path_to_meta "git_revision.m"], filesep));

catch
    % nothing to do
    warning("Not inside a git repository.")
end

end
