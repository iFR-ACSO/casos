% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function f = to_function(obj,variables,opts)
% Return casadi.Function object using SX.

if (nargin < 2) || isstruct(variables)
    % default variables
    if (nargin > 1) 
        % options given
        opts = variables;
    else
        % default options
        opts = struct;
    end
    variables = tuple2cell(obj.indeterminates);

elseif (nargin < 3)
    % default options
    opts = struct;
end

% create symbolic variables
in = cellfun(@to_symbolic, variables, 'UniformOutput', false);

% substitute symbolic variables
p = subs(obj,[variables{:}],casos.package.polynomial(vertcat(in{:})));

f = casadi.Function('f',in,{casadi.SX(p)},opts);

end

function var = to_symbolic(x)
% Convert indeterminate variable(s) into SX symbol.

names = str(x);

% use first character as name
var = casadi.SX.sym(names{1}(1), length(x), 1);

end
