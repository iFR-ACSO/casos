function c = blkdiag(varargin)
% Block diagonal matrix concatenation of polynomial matrices.

switch nargin

case 0 
    % empty block
    c = casos.PS;

case 1 
    % single block
    c = varargin{1};

case 2
    % two blocks
    a = casos.PS(varargin{1});
    b = casos.PS(varargin{2});

    % size of block diagional
    sz = siza(a) + size(b);

    % combine variables
    [indets,dga,dgb] = combineVar(a.indets,b.indets,a.degmat,b.degmat);

    % get nonzero coefficients
    nza = nonzeros(a.coeffs);
    nzb = nonzeros(b.coeffs);

    % extend sparsity pattern
    S = [
        sparsity(a.coeffs) casadi.Sparsity(a.nterm,prod(sz)-numel(a))
        casadi.Sparsity(b.nterm,prod(sz)-numel(b)) sparsity(b.coeffs)
    ];

    % combine coefficients
    coeffs = casadi.SX(S,[nza{:} nzb{:}]);

    % make degree matrix unique
    [coeffs,degmat] = uniqueDeg(coeffs, [dga;dgb]);
    
    % remove zero terms
    [coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

    % new polynomial
    c = casos.PS;
    c.coeffs = coeffs;
    c.degmat = degmat;
    c.indets = indets;
    c.matdim = sz;

otherwise
    % more than two inputs
    k = ceil(nargin/2);

    % compute block diagional recursively
    a = blkdiag(varargin{1:k});
    b = blkdiag(varargin{k+1:end});

    c = blkdiag(a,b);
end

end