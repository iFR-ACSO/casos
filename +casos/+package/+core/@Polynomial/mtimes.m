% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = mtimes(a,b)
% Matrix multiplication of two polynomials.

if isscalar(a) || isscalar(b)
    % fall back to scalar multiplication
    c = times(a,b);
    return

elseif ~check_sz_matrix(a,b)
    % dimensions are compatible if size(a,2) == size(b,1)
    throw(casos.package.core.IncompatibleSizesError.matrix(a,b));

elseif isempty(a) || isempty(b)
    % empty multiplication
    c = a.zeros(size(a,1),size(b,2));
    return
end

% else
c = a.new_poly;

% compute coefficient matrix and sparsity
[S,c.coeffs] = mtimes_internal(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
