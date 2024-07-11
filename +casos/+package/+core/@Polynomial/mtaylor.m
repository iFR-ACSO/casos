function c = mtaylor(a,x,b,deg)
% Symbolic Taylor expansion.

assert(nargin == 4, 'Incorrect number of inputs (got: %d, expected 4).',nargin)

assert(is_symbolic(x), 'Second argument must be purely symbolic.')

% check dimensions
if isscalar(b)
    % use same expression
    b = repmat(b,size(x));

else
    assert(numel(x) == numel(b),'Second and third argument have incompatible sizes.')
end

% else:
c = a.new_poly;

% project to basis
[X,S] = poly2basis(x);
[B] = poly2basis(b,S);

% linearize coefficients
coeffs = mtaylor(a.coeffs,X,B,deg);

% remove zero terms
[S,c.coeffs] = coeff_update(a.get_sparsity,coeffs);

% set sparsity
c = set_sparsity(c,S);

end

