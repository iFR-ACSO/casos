% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = times(a,b)
% Element-wise multiplication of two polynomial matrices.

% input dimensions
sza = size(a);
szb = size(b);

% dimensions are compatible if equal or one summand is row/column
[tf,sz] = check_sz_basic(a,b);

if (~tf)
    throw(casos.package.core.IncompatibleSizesError.basic(a,b));
end

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % element-wise multiplication with empty polynomial is empty
    c = a.empty(sz);
    return
end
% TODO: handle or escape for other simple cases, e.g., scalar, constant
% matrix, single term etc.?

% else
c = a.new_poly;

% reshape to output dimensions
[S1,cfa] = coeff_repmat(a.get_sparsity,a.coeffs,sz./sza);
[S2,cfb] = coeff_repmat(b.get_sparsity,b.coeffs,sz./szb);

% multiply coefficient matrices
[S,c.coeffs] = coeff_times(S1,S2,cfa,cfb);

% set sparsity
c = set_sparsity(c,S);

end
