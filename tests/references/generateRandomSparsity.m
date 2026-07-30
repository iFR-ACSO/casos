% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function sp = generateRandomSparsity(sz,variables,degmax,density,monomial)
% Generate a random sparsity pattern.

select_variables = (randn(length(variables),1) > 0);

% select random variables
x = variables(select_variables);

% number of variables
n = length(x);

% create subindices (0-index)
I = 0:sz(1)-1;
J = 0:sz(2)-1;

% permute subindices
[I,J] = ndgrid(I,J);

% create linear indices
idx = 1:numel(I);

if (nargin > 3 && density < 1)
    % generate matrix sparsity
    idx(rand(size(idx)) > density) = 0;
end

if (nargin > 4 && monomial)
    % create matrix of monomials
    D = arrayfun(@(~) randi(degmax,[length(idx) 1]), 1:n, 'UniformOutput', false);
    % select all nonzero indices
    select_entries = (idx(:) ~= 0);

else
    % create degree vector
    D = cell(1,n);
    D(:) = {0:degmax};
    
    % permute linear indices and degrees
    [idx,D{:}] = ndgrid(idx,D{:});
    
    % select random indices and degrees
    select_entries = (idx(:) ~= 0 & randn(numel(idx),1) > 0);
end

idx = squeeze(idx(select_entries));
i = squeeze(I(idx));
j = squeeze(J(idx));
degree_vectors = cellfun(@(d) squeeze(d(select_entries)), D, 'UniformOutput', false);

if isempty(degree_vectors)
    % empty degrees
    degrees = zeros(length(i),length(x));
else
    % concatenate
    degrees = reshape([degree_vectors{:}],length(i),length(x));
end

% create sparsity pattern
sp = casos.Sparsity.tuplet(sz(1),sz(2),i,j,x,degrees);

end
