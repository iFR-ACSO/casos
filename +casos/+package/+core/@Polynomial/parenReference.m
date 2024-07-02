function varargout = parenReference(obj,indexOp)
% Overwriting matlab.mixin.indexing.RedefinesParen.parenReference

% perform parenthesis reference
idx = indexOp(1);

if length(indexOp) > 1
    % handle getters on referenced polynomial
    [varargout{1:nargout}] = parenReference@casos.package.core.GenericPolynomial(obj,indexOp);
    return

elseif isequal(idx.Indices,{':'})
    % handle vectorization
    varargout = {reshape(obj,numel(obj),1)};
    return

elseif isequal(idx.Indices,{':' ':'})
    % nothing to do
    varargout = {obj};
    return
end

% select referenced elements
I = reshape(1:numel(obj),size(obj));
ii = I.(idx);

if ~isempty(ii)
    % new polynomial
    p = obj.new_poly;
    % reference coefficients
    [S,p.coeffs] = coeff_subsref(obj.get_sparsity,obj.coeffs,ii,size(ii));

    % set sparsity
    p = set_sparsity(p,S);
else
    % empty reference
    p = obj.zeros(size(ii));
end

% return
varargout = {p};

end
