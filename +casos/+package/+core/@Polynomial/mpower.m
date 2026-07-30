% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = mpower(a,n)
% Power of polynomial matrix.

if ~check_sz_square(a) || ~isscalar(n)
    % dimensions are compatible if a is square and n is scalar
    throw(casos.package.core.IncompatibleSizesError.mpower(a,n));

elseif isscalar(a)
    % fall back to scalar power.
    b = power(a,n);
    return

elseif n == 0
    % identity matrix
    b = a.eye(length(a));
    return

elseif n == 1
    % no change
    b = a;
    return
end

b = a.new_poly;

% Call mpower recursively
[S,b.coeffs] = mpower_internal(a.get_sparsity,a.coeffs,n);

% Note: A for loop on mtimes requires n calls to mtimes. A recursion
% only requires roughly log2(n) calls to mtimes.

b = set_sparsity(b,S);

end
