function obj = parenAssign(obj,indexOp,varargin)
% Overwriting matlab.mixin.indexing.RedefinesParen.parenAssign

idx = indexOp(1);

if length(indexOp) > 1
    % perform parantheses reference
    p = obj.(idx);
    % forward assign
    [p.(indexOp(2:end))] = varargin{:};
    % re-assign modified element
    obj.(idx) = p;

    return
end

% else
assert(length(varargin) == 1, 'Too many arguments on right-hand side.')

% assignment
q = casos.PS(varargin{:});

nta = obj.nterm;
ntb = q.nterm;

% combine variables
[indets,dga,dgb] = combineVar(obj.indets,q.indets,obj.degmat,q.degmat);

% select referenced elements
I = sparse(size(obj,1),size(obj,2));
I.(idx) = 1;

% size of reference
sz = size(I.(idx));

% reshape to referenced size
cfb = reshape(repmat(q.coeffs,sz./size(q)),ntb,prod(sz));

% extend coefficient matrix
S = sparsity(casadi.SX(repmat(I(:)',ntb,1)));
% assign rhs coefficients
coeffs = [obj.coeffs; casadi.SX(S,cfb)];

% remove selected lhs coefficients
coeffs(1:nta,find(I)) = 0;
% coeffs(nta+(1:ntb),find(I)) = cfb;

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, [dga;dgb]);

% remove zero terms
[obj.coeffs,obj.degmat,obj.indets] = removeZero(coeffs,degmat,indets);

end
