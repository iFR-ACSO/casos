function obj = new_sym(obj,dstr,varargin)
% Create a new symbolic polynomial (helper function).

if nargin < 3
    % constant scalar
    S = casos.Sparsity.scalar;

elseif ~isa(varargin{1},'casos.Sparsity')
    % zero-degree symbol
    S = casos.Sparsity(varargin{:});

elseif nargin > 3
    % repeat scalar pattern
    assert(isscalar(varargin{1}),'Sparsity pattern must be scalar.')

    S = casos.dense(varargin{:});

else
    % sparsity pattern
    S = varargin{1};
end

% assign sparsity
obj = set_sparsity(obj,S);
% symbolic coefficient matrix
obj.coeffs = sym_coeff(obj,dstr,coeff_sparsity(S));

end
