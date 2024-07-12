function p = cat(dim,varargin)
% Concatenate polynomials along specified dimension.
%
% Overwriting matlab.mixin.indexing.RedefinesParen.cat

if nargin ~= 3
    % handle non-binary concatenation
    p = cat@casos.package.core.PolynomialInterface(dim,varargin{:});
    return
end

% else: two polynomials given
p1 = varargin{1};
p2 = varargin{2};

switch (dim)
    case 0
        % block diagonal
        ul = p1.zeros(size(p1,1),size(p2,2));
        rl = p1.zeros(size(p2,1),size(p1,2));

        % concatenate all blocks
        p = blockcat(p1,ul,rl,p2);

    otherwise
        % generic concatenation
        p = p1.new_poly;

        % concatenate coefficient matrices
        [S,p.coeffs] = coeff_cat(p1.get_sparsity,p2.get_sparsity,p1.coeffs,p2.coeffs,dim);
        % set sparsity
        p = set_sparsity(p,S);
end

end
