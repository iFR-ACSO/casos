% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = power(a,n)
% Element-wise power of polynomial sparsity pattern.

a = casos.Sparsity(a);

% sparsity pattern and exponent must be of same size
if ~check_sz_equal(a,n)
    throw(casos.package.core.IncompatibleSizesError.other(a,n));

elseif all(n==0)
    % power of zero is dense
    b = casos.Sparsity.dense(size(a));
    return
end

% power of coefficient matrix
b = coeff_power(a,a.coeffs,n);

end
