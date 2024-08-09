function c = linearize(a,x,b)
% Symbolic linearization of polynomial expression.

assert(nargin == 3, 'Incorrect number of inputs (got: %d, expected 3).',nargin)

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
coeffs = linearize(a.coeffs,X,B);

% remove zero terms
[S,c.coeffs] = coeff_update(a.get_sparsity,coeffs);

% set sparsity
c = set_sparsity(c,S);

end
