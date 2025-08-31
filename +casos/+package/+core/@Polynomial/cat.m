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

% generic concatenation
p = p1.new_poly;

% concatenate coefficient matrices
[S,p.coeffs] = coeff_cat(p1.get_sparsity,p2.get_sparsity,p1.coeffs,p2.coeffs,dim);
% set sparsity
p = set_sparsity(p,S);

end
