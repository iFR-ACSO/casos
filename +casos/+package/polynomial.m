function p = polynomial(a,varargin)
% Convert algebraic input to polynomial.

if isa(a,'casos.package.core.GenericPolynomial')
    % nothing to do
    assert(nargin < 2,'Too many input arguments.')

    p = a;
    return

elseif isa(a,'casos.Indeterminates')
    % indeterminate variables
    p = casos.PD(a);
    return

elseif isa(a,'casos.Sparsity')
    % first input is sparsity pattern
    assert(nargin > 1,'Not enough input arguments.')
    assert(nargin < 3,'Too many input arguments.')

    S = {a};
    a = varargin{1};

else
    % single input
    S = {};
end

% choose suitable polynomial type
switch class(a)
    case {'double' 'casadi.DM'}
        % double polynomial
        p = casos.PD(S{:},a);

    case 'casadi.SX'
        % symbolic polynomial
        p = casos.PS(S{:},a);

    case 'casadi.MX'
        % symbolic matrix polynomial
        error('Not supported.')

    otherwise
        error('No conversion from %s to any polynomial class.', class(a))
end

end
