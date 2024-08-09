function S = cat(dim,varargin)
% Concatenate polynomial sparsity patterns along specified dimension.
%
% Overwriting matlab.mixin.indexing.RedefinesParen.cat

if nargin ~= 3
    % handle non-binary concatenation
    p = cat@casos.package.core.PolynomialInterface(dim,varargin{:});
    return
end

% two patterns given
S1 = casos.Sparsity(varargin{1});
S2 = casos.Sparsity(varargin{2});

% concatenate coefficient matrices
S = coeff_cat(S1,S2,S1.coeffs,S2.coeffs,dim);

end
