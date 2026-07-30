% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = dot(a,b)
% Scalar dot product of two polynomials.

if ~check_sz_equal(a,b)
    % dimension mismatch
    throw(casos.package.core.IncompatibleSizesError.other(a,b));
end

% expand coefficient matrices to match
[cf1,cf2] = coeff_expand(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% return dot product of coefficients
c = dot(cf1,cf2);

end
