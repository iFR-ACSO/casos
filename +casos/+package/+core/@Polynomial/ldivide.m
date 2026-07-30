% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = ldivide(a,p)
% Element-wise left array division with constant matrix.

assert(is_zerodegree(a),'Only division by constant or symbolic matrix possible.')

% input dimensions
sza = size(a);
szp = size(p);

% dimensions are compatible if equal or one summand is row/column
[tf,sz] = check_sz_basic(a,p);

if (~tf)
    throw(casos.package.core.IncompatibleSizesError.basic(a,p));
end

% handle simple case(s) for speed up
if isempty(a) || isempty(p)
    % element-wise multiplication with empty polynomial is empty
    c = p.empty(sz);
    return
end

% else
c = p.new_poly;

% reshape to output dimensions
[Sa,cfa] = coeff_repmat(a.get_sparsity,a.coeffs,sz./sza);
[Sp,cfp] = coeff_repmat(p.get_sparsity,p.coeffs,sz./szp);

% element-wise division of coefficients
coeffs = coeff_repterms(Sa,cfa,p.nterm) .\ cfp;

% update sparsity pattern
[S,c.coeffs] = coeff_update(Sp,coeffs,sz);

% set sparsity
c = set_sparsity(c,S);

end
