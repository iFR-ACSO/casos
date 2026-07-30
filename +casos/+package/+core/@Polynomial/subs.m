% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = subs(a,x,b)
% Polynomial substitution of indeterminates x in a by expression b.

assert(is_indet(x),'Second argument must be vector of indeterminate variables.')

% check dimensions
if isscalar(b)
    % replace all indeterminates by same expression
    b = repmat(b,size(x));

elseif (size(x,1) == size(b,1) && iscolumn(x) && size(b,2) > 1)
    % repeated substitution -- not supported
    error('Repeated substitution not supported, use to_function(a) instead.')

elseif ~isequal(numel(x), numel(b))
    % number of elements to be substituted is incompatible
    throw(casos.package.core.IncompatibleSizesError.substitute(x,b));
end

c = a.new_poly;

% substitute indeterminates in coefficient matrix
[S,c.coeffs] = coeff_subs(a.get_sparsity,a.coeffs,x,b.get_sparsity,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
