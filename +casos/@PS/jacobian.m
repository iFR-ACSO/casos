function g = jacobian(f,x)
% Compute Jacobian matrix of vector polynomial expression.

assert(is_indet(x), 'Second argument must be vector of indeterminates.')

% else:
f = casos.PS(f);

% number of terms in f
nt = f.nterm;

% number of elements in f
ne = f.numel;

% number of indeterminates
nx = length(x);

% find position of x(i) in the degree matrix of f
[tf,I] = ismember(x.indets,f.indets);

% jacobian(sum_a c_a*x^a1...x^an, x) = [
%       sum_a a1*c_a*x^(a1-1)...x^an
%       ...
%       sum_a an*c_a*x^a1*...*x^(an-1)
% ]
g = casos.PS;

% extend coefficient + degree matrices to ne*nx elements
coeffs = kron(x.degmat*x.coeffs,f.coeffs);
degmat = kron(ones(nx,1),f.degmat);
% indices to degrees per indeterminates
I2 = kron(I(:,tf),ones(nt,1));

% multiply by degree associated with ai
scale = zeros(nt*nx,1);
% enumerate extended coefficients
idx = reshape(1:(nt*nx),nt,nx);
% remove indeterminates f is constant in
idx(:,~tf) = [];
% index into degree matrix
ind = sub2ind(size(degmat),idx,I2);

% get multipliers
scale(idx) = degmat(ind);
% decrease corresponding degrees by 1
degmat(ind) = degmat(ind)-1;

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(scale.*coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,f.indets);

% new polynomial
g.coeffs = coeffs;
g.degmat = degmat;
g.indets = indets;
g.matdim = [ne nx];

end
