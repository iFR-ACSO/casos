function obj = parenAssign(obj,indexOp,varargin)
% Implementing matlab.mixin.indexing.RedefinesParen.parenAssign

idx = indexOp(1);

if length(indexOp) > 1
    % perform parentheses reference
    v = obj.(idx);
    % forward assign
    [p.(indexOp(2:end))] = varargin{:};
    % re-assign modified element
    obj.(idx) = v;

    return
end

% else
assert(length(varargin) == 1, 'Too many arguments on right-hand side.')

% assignment
args = casos.Indeterminates(varargin{:});

% get LHS variables
vars = obj.variables;

% assign RHS variables
vars.(idx) = args.variables;

% return
obj = casos.Indeterminates(vars{:});

end
