% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function spy(obj)
% Visualize polynomial sparsity pattern.

[degree,L] = get_degree(obj);
% repeat degree vector to match logical matrix
degree = repmat(degree,numel(obj),1);

% degrees of each term
degrees = spalloc(size(L,1),size(L,2),nnz(L));
degrees(L) = degree(L);

% element-wise max degree
maxdeg = max(degrees,[],2);

% detect nonzero elements
Lnz = any(L,2);

% prepare output
out = cell(size(obj));

% set per-element output (sparse zero = ., constant = *, degree = k)
out(~Lnz) = {'.'};
out(Lnz & maxdeg == 0) = {'*'};
out(Lnz & maxdeg > 9) = {'x'};
out(maxdeg > 0 & maxdeg < 10) = compose('%d',full(maxdeg(maxdeg > 0 & maxdeg < 10)));

% print output to command line
disp(cell2mat(out))
disp(' ')

end
