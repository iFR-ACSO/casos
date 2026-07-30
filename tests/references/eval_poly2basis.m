% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function res = eval_poly2basis(arg1,arg2)
% Evaluate poly2basis on first argument with basis of second argument.

assert(isequal(size(arg1),size(arg2)), 'Size mismatch.')

if isempty(arg1) || isempty(arg2)
    % handle empty polynomial or basis
    res = zeros(0,1);
    return
end

% else
N = prod(size(arg1)); %#ok<PSIZE>

% evaluate poly2basis element-wise
results = arrayfun(@(i) full(poly2basis(arg1(i),monomials(arg2(i)))), 1:N, "UniformOutput",false);

% concatenate results to column vector
res = vertcat(results{:});

end
