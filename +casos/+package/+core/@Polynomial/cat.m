% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = cat(dim,varargin)
% Concatenate polynomials along specified dimension.
%
% Overwriting matlab.mixin.indexing.RedefinesParen.cat

if nargin ~= 3
    % handle non-binary concatenation
    c = cat@casos.package.core.PolynomialInterface(dim,varargin{:});
    return
end

% else: two polynomials given
a = varargin{1};
b = varargin{2};

switch (dim)
    case 0
        % block diagonal
        ul = a.zeros(size(a,1),size(b,2));
        rl = a.zeros(size(b,1),size(a,2));

        % concatenate all blocks
        c = blockcat(a,ul,rl,b);

    otherwise
        if (isempty(a) && isempty(b))
            % concatenation of empty polynomials
            sza = size(a);
            szb = size(b);

            if (sza(3-dim) == 0 && szb(3-dim) == 0)
                % concatenate along non-empty dimension
                sz = sza + szb;
            else
                % result is empty along dimension
                sz(dim) = min(sza(dim), szb(dim));
                sz(3-dim) = max(sza(3-dim), szb(3-dim));
            end
            c = a.empty(sz);
            return

        elseif isempty(a)
            % concatenation of b with empty polynomial
            c = b;
            return
        
        elseif isempty(b)
            % concatenation of a with empty polynomial
            c = a;
            return

        elseif (size(a,3-dim) ~= size(b,3-dim))
            % dimensions are consistent if size(a,~dim) == size(b,~dim)
            throw(casos.package.core.IncompatibleSizesError.concat(a,b));
        end

        % generic concatenation
        c = a.new_poly;

        % concatenate coefficient matrices
        [S,c.coeffs] = coeff_cat(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs,dim);
        % set sparsity
        c = set_sparsity(c,S);
end

end
