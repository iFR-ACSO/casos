function varargout = parenReference(obj,indexOp)
% Implementing matlab.mixin.indexing.RedefinesParen.parenReference

idx = indexOp(1);

% perform reference
vars = obj.variables.(idx);

% new indeterminates
out = casos.Indeterminates;
out.variables = vars;

if length(indexOp) > 1
    % forward reference
    [varargout{1:nargout}] = out.(indexOp(2:end));

else
    varargout = {out};
end

end
