% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = cat(dim,varargin)
% Concatenate polynomial sparsity patterns along specified dimension.
%
% Overwriting matlab.mixin.indexing.RedefinesParen.cat

if nargin ~= 3
    % handle non-binary concatenation
    c = cat@casos.package.core.PolynomialInterface(dim,varargin{:});
    return
end

% two patterns given
a = casos.Sparsity(varargin{1});
b = casos.Sparsity(varargin{2});

[tf,sz] = check_sz_concat(dim,a,b);

if (~tf)
    % dimensions are consistent if size(a,~dim) == size(b,~dim)
    throw(casos.package.core.IncompatibleSizesError.concat(a,b));

elseif (isempty(a) && isempty(b))
    % concatenation of empty polynomials
    c = casos.Sparsity(sz);
    return

elseif isempty(a)
    % concatenation of b with empty polynomial
    c = b;
    return

elseif isempty(b)
    % concatenation of a with empty polynomial
    c = a;
    return

else
end

% concatenate coefficient matrices
c = coeff_cat(a,b,a.coeffs,b.coeffs,dim);

end
