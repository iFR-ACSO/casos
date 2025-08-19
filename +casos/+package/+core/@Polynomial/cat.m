function p = cat(dim,varargin)
% Concatenate polynomials along specified dimension.
%
% Overwriting matlab.mixin.indexing.RedefinesParen.cat

if nargin ~= 3
    % handle non-binary concatenation
    p = cat@casos.package.core.PolynomialInterface(dim,varargin{:});
    return
end

% else: two vectors given
p1 = varargin{1};
p2 = varargin{2};

assert(is_operator(p1) == is_operator(p2), 'Must not mix polynomials and operators.')

switch (dim)
    case 0
        % block diagonal
        if is_operator(p1)
            % fill with zero operators
            ur = p1.zero_operator(size(p1,1),size(p2,2));
            ll = p1.zero_operator(size(p2,1),size(p1,2));
        else
            % fill with zero matrices
            ur = p1.zeros(size(p1,1),size(p2,2));
            ll = p1.zeros(size(p2,1),size(p1,2));
        end

        % concatenate all blocks
        p = blockcat(p1,ur,ll,p2);

    otherwise
        % generic concatenation
        p = p1.new_poly;

        % concatenate coefficient matrices
        [S,p.coeffs] = coeff_cat(p1.get_sparsity,p2.get_sparsity,p1.coeffs,p2.coeffs,dim);
        % set sparsity
        p = set_sparsity(p,S);
end

end
