function obj = parenDelete(obj,indexOp)
% Overwriting matlab.mixin.indexing.RedefinesParen.parenDelete

% MATLAB uses different indexing type for delete
i = indexOp(1).Indices;

if length(indexOp) > 1
    % perform parantheses reference
    p = obj(i{:});
    % forward delete
    p.(indexOp(2:end)) = [];
    % re-assign modified element
    obj(i{:}) = p;

    return
end

% else
I = sparse(size(obj,1),size(obj,2));
I(i{:}) = 1;

% removed referenced coefficients
coeffs = obj.coeffs(:,find(~I));

% remove coefficients, degrees, and/or indeterminates 
% that do not appear in the referenced polynomial
[obj.coeffs,obj.degmat,obj.indets] = removeZero(coeffs,obj.degmat,obj.indets);
% resize
I(i{:}) = [];
obj.matdim = size(I);

end
