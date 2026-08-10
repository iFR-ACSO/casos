% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (Sealed) CasosMeta
% Meta information for CaΣoS.
%
% - version: Returns version in the format 'x.y.z(+)'.
% - git_revision: Returns git revision ID.
%

methods
    function obj = CasosMeta()
        error('Forbidden.')
    end
end

methods (Static)
    v = version;
    v = git_revision;
end

end
