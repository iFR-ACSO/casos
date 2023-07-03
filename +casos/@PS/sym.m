function p = sym(dstr,w,sz,type)
% Create symbolic polynomial.

if nargin < 4
    % default type: polynomial
    type = '';
end
if nargin < 3
    % default size: scalar
    sz = [1 1];
elseif isscalar(sz)
    sz = [sz sz];
elseif ~isrow(sz) || ~numel(sz) == 2
    error('Third input must be scalar or 1x2 vector of dimensions.')
end
if nargin < 2
    % default monomials: 1
    w = casos.PS(1);
else
    w = casos.PS(w);
    assert(is_monom(w),'Second input must be a vector of monomials.')
end
if nargin < 1
    error('Undefined inputs.');
end

p = casos.PS;

% number of monomial terms
nt = w.nterm;
% number of elements requested
ne = prod(sz);

switch type
    case 'gram'
        % create using Gram matrix form
        Q = casadi.SX.sym(dstr,nt^2,ne);  % Gram matrices, TODO: symmetric?
        D = kron(w.degmat,ones(nt,1)) + kron(ones(nt,1),w.degmat);
        % make degree matrix unique
        [p.degmat,ia,ic] = unique(D,'rows','stable');
        % sum repeated coefficients
        summat = sparse(ic,1:(nt^2),1,length(ia),nt^2);
        p.coeffs = summat*Q;

    otherwise
        % create using coefficient vector form
        p.coeffs = casadi.SX.sym(dstr,nt,ne);
        p.degmat = w.degmat; % TODO: match ordering in w
end

% set indeterminates + dimensions
p.indets = w.indets;
p.matdim = sz;

end