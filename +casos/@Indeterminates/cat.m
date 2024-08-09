function v = cat(dim,varargin)
% Concatenate indeterminate variables.
%
% Note: Indeterminate variables are always concatenated horizontally.

% parse input arguments
[tf,vars] = cellfun(@parse, varargin, 'UniformOutput',false);

if ~all([tf{:}])
    % concatenate algebraic objects
    v = cat@casos.package.core.AlgebraicObject(dim,varargin{:});
    return
end

% else:
vars = [vars{:}];

% concatenate variables
v = casos.Indeterminates(vars{:});

end

function [tf,vars] = parse(arg)
% Parse input arguments.

% check for indeterminate variables
tf = (isa(arg,'casos.package.core.AlgebraicObject') && is_indet(arg));

if ~tf
    % nothing to do
    vars = {};
else
    % get variables
    arg = casos.Indeterminates(arg);

    vars = arg.variables;
end

end
