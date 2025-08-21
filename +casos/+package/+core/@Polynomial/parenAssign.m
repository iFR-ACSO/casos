function obj = parenAssign(obj,indexOp,varargin)
% Overwriting matlab.mixin.indexing.RedefinesParen.parenAssign

idx = indexOp(1);

if length(indexOp) > 1
    % handle forwarded assignments to subreference
    obj = parenAssign@casos.package.core.GenericPolynomial(obj,indexOp,varargin{:});
    return
end

% else
assert(nargin < 4, 'Too many arguments on right-hand side.')

% assignment
q = obj.new_poly(varargin{:});

assert(~is_operator(obj) && ~is_operator(q), 'Subscripted assignment of operators not supported.')

% select referenced elements
I = sparse(size(obj,1),size(obj,2));
I.(idx) = 1;

% reference
lhs = I.(idx);

% dimensions are compatible if equal or right side is row/column
if ~check_sz_assign(lhs,q)
    throw(casos.package.core.IncompatibleSizesError.other(lhs,q))
end

% reshape coefficients to referenced size
[S2,cfb] = coeff_repmat(get_sparsity(q),q.coeffs,size(lhs)./size(q));

% assign coefficients
[S,obj.coeffs] = coeff_subsasgn(obj.get_sparsity,S2,obj.coeffs,cfb,find(I));

% set sparsity
obj = set_sparsity(obj,S);

end
