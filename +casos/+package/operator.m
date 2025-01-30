function op = operator(a,varargin)
% Convert algebraic input to linear operator.

if isa(a,'casos.package.core.OperatorSparsity')
    % first input is operator sparsity pattern
    assert(nargin > 1,'Not enough input arguments.')
    assert(nargin < 3,'Too many input arguments.')

    S = {a};
    a = varargin{1};
    varargin = {};

else
    % single input
    S = {};
end

% choose suitable operator type
switch class(a)
    case {'double' 'casadi.DM'}
        % double polynomial operator
        op = casos.PDOperator(S{:},a,varargin{:});

    case 'casadi.SX'
        % symbolic polynomial operator
        op = casos.PSOperator(S{:},a,varargin{:});

    case 'casadi.MX'
        % symbolic matrix polynomial operator
        error('Not supported.')

    otherwise
        error('No conversion from %s to any polynomial class.', class(a))
end

end
