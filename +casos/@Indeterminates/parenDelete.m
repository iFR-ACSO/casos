function obj = parenDelete(obj,indexOp)
% Implementing matlab.mixin.indexing.RedefinesParen.parenDelete

% MATLAB uses different indexing type for delete
i = indexOp(1).Indices;

if length(indexOp) > 1
    % perform parantheses reference
    v = obj(i{:});
    % forward delete
    v.(indexOp(2:end)) = [];
    % re-assign modified element
    obj(i{:}) = v;

    return
end

% else
vars = obj.variables;

% remove referenced variables
vars(i{:}) = [];

% return
obj.variables = vars;

end
