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

    % sizes
    sza = size(a);
    szb = size(b);

    % terms
    nta = a.nterm;
    ntb = b.nterm;

    % size of block diagonal
    sz = sza + szb;

    % combine variables
    [indets,dga,dgb] = combineVar(a.indets,b.indets,a.degmat,b.degmat);

    % reshape to access columns
    cfa = prepareMatrixOp(a.coeffs,sza,1);
    cfb = prepareMatrixOp(b.coeffs,szb,1);

    % add zeros to each column
    cfA = [cfa; sparse(szb(1),size(cfa,2))];
    cfB = [sparse(sza(1),size(cfb,2)); cfb];

    % combine coefficients
    coeffs = [
        finishMatrixOp(cfA,[sz(1) sza(2)],1) sparse(nta,sz(1)*szb(2))
        sparse(ntb,sz(1)*sza(2)) finishMatrixOp(cfB,[sz(1) szb(2)],1)
    ];

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