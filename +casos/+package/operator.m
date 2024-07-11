function op = operator(a,varargin)
% Convert algebraic input to linear operator.

% choose suitable operator type
switch class(a)
    case {'double' 'casadi.DM'}
        % double polynomial
        op = casos.PDOperator(a,varargin{:});

    case 'casadi.SX'
        % symbolic polynomial
        op = casos.PSOperator(a,varargin{:});

    case 'casadi.MX'
        % symbolic matrix polynomial
        error('Not supported.')

    otherwise
        error('No conversion from %s to any polynomial class.', class(a))
end

end
