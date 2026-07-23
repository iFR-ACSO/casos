% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function S = basis(p,I)
% Return basis of multidimensional polynomial.
%
% Note: This function is made obsolete by the use of polynomial sparsity
% and only kept for use with indexing argument.

if nargin < 2
    % fall back to sparsity
    warning('Deprecated: Use sparsity() instead.')
    S = sparsity(p);

elseif (numel(p) ~= length(I))
    % dimension mismatch
    error('Dimension mismatch.')

elseif isempty(p) || all(~I)
    % return empty basis
    S = casos.Sparsity(0,1);

else
    % return sparsity of subindex
    S = sub(p.get_sparsity,find(I),casos.Sparsity(nnz(I),1));
end

end
