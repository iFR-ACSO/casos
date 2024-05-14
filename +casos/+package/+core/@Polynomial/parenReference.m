function varargout = parenReference(obj,indexOp)
% Overwriting matlab.mixin.indexing.RedefinesParen.parenReference

% perform parenthesis reference
idx = indexOp(1);

% select referenced elements
I = logical(sparse(size(obj,1),size(obj,2)));
I.(idx) = true;

if length(indexOp) > 1
    % handle getters on referenced polynomial
    [varargout{1:nargout}] = parenReference@casos.package.core.GenericPolynomial(obj,indexOp);
    return
end

% new polynomial
p = obj.new_poly;

if nnz(I) > 0
    % reference coefficients
    [S,p.coeffs] = coeff_subsref(obj.get_sparsity,obj.coeffs,find(I),size(I.(idx)));

    % set sparsity
    p = set_sparsity(p,S);
end

% return
varargout = {p};

end
