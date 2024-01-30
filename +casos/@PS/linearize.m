function c = linearize(a,x,b)
% Symbolic linearization of polynomial expression.

assert(nargin == 3, 'Incorrect number of inputs (got: %d, expected 3).',nargin)

a = casos.PS(a);
x = casos.PS(x);
b = casos.PS(b);

assert(is_symbolic(x), 'Second argument must be purely symbolic.')

% check dimensions
if isscalar(b)
    % use same expression
    b = repmat(b,size(x));

else
    assert(numel(x) == numel(b),'Second and third argument have incompatible sizes.')
end

% else:
c = casos.PS;

% project to basis
[Q,z] = poly2basis(x);
[P] = poly2basis(b,z);

% linearize coefficients
coeffs = linearize(a.coeffs,Q,P);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,a.degmat,a.indets);

% new polynomial
c.coeffs = coeffs;
c.degmat = degmat;
c.indets = indets;
c.matdim = a.matdim;

end
