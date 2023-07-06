function varargout = parenReference(obj,indexOp)
% Overwriting matlab.mixin.indexing.RedefinesParen.parenReference

% perform parenthesis reference
idx = indexOp(1);

% select referenced elements
I = sparse(size(obj,1),size(obj,2));
I.(idx) = 1;

% reference coefficients
coeffs = obj.coeffs(:,find(I));

% new polynomial
p = casos.PS;
% remove coefficients, degrees, and/or indeterminates 
% that do not appear in the referenced polynomial
[p.coeffs,p.degmat,p.indets] = removeZero(coeffs,obj.degmat,obj.indets);
% resize
p.matdim = size(I.(idx));


if length(indexOp) > 1
    % forward reference
    [varargout{1:nargout}] = p.(indexOp(2:end));

else
    varargout = {p};
end

end
