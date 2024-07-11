function [S,coeffs] = coeff_nabla(S,coeffs,x)
% Compute coefficient matrix for polynomial nabla operator.

x = casos.Indeterminates(x);

% number of terms in f
nt = S.nterm;

% number of elements in f
ne = S.numel;

% number of indeterminates
nx = length(x);

% find position of x(i) in the degree matrix of f
[tf,I] = ismember(x,S.indets);

% adjust for order of variables in x
x_coeffs = casadi.Sparsity.diag(nx);

% extend coefficient + degree matrices to ne*nx elements
coeffs = kron(x_coeffs,coeffs);
degmat = kron(ones(nx,1),S.degmat);
% indices to degrees per indeterminates
I2 = kron(I(:,tf),ones(nt,1));

% multiply by degree associated with ai
scale = zeros(nt*nx,1);
% enumerate extended coefficients
idx = reshape(1:(nt*nx),nt,nx);
% remove indeterminates S is constant in
idx(:,~tf) = [];
% index into degree matrix
ind = sub2ind(size(degmat),idx,I2);

% get multipliers
scale(idx) = degmat(ind);
% decrease corresponding degrees by 1
degmat(ind) = degmat(ind)-1;

if ~isa(coeffs,'casadi.Sparsity')
    % scale coefficients
    coeffs = scale.*coeffs;
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,S.indets);

% return sparsity pattern
S.degmat = degmat;
S.indets = indets;
S.matdim = [ne nx];
% store coefficients
S = set_coefficients(S,coeffs);

end
