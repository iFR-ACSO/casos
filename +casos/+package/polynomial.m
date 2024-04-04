function p = polynomial(a)
% Convert algebraic input to polynomial.

if isa(a,'casos.package.core.AlgebraicObject')
    % nothing to do
    p = a;
    return
end

% else
switch class(a)
    case {'double' 'casadi.DM'}
        % double polynomial
        p = casos.PS(a);

    case 'casadi.SX'
        % symbolic polynomial
        p = casos.PS(a);

    case 'casadi.MX'
        % symbolic matrix polynomial
        error('Not supported.')

    otherwise
        error('No conversion from %s to any polynomial class.', class(a))
end

end
