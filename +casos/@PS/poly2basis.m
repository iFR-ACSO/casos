function [q,z,e] = poly2basis(p,z,I,legacy)
% Project polynnomial p onto the span of monomials contained in z.

p = casos.PS(p);

if nargin < 2
    % syntax poly2basis(p)
    z = basis(p);
    I = true(size(p));
    legacy = true;
elseif islogical(z)
    % syntax poly2basis(p,I)
    I = z;
    z = basis(p,I);
    legacy = true;
elseif nargin < 3
    % syntax poly2basis(p,z,I)
    I = true(size(p));
    legacy = true;
elseif isempty(z) || isempty(I)
    % syntax poly2basis(...,legacy)
    if isempty(I), I = true(size(p)); end
    if isempty(z), z = basis(p,I); end
end

% combine variables
[~,dgp,dgz] = combineVar(p.indets,z.indets,p.degmat,z.degmat);

% find degrees of z in p
[tf,ii] = ismember(dgp,dgz,'rows');
idx = 1:p.nterm;

if iscolumn(z) && all(I) && legacy
    % legacy code: project onto monomial vector
    assert(is_monom(z),'Second argument must be vector of monomials.')

    % select coefficients
    q = casadi.SX.zeros(numel(z),numel(p));

    if all(~tf)
        % zero projection
        e = p;
        return
    end

    % else
    q(ii(tf),:) = p.coeffs(idx(tf),:);
    % reorder to match monomials
    q = z.coeffs * q;

else
    % project onto basis matrix
    [lZ,lp] = size(z);

    if isscalar(p)
        cfp = repmat(p.coeffs,1,lp);
    elseif lp > 0
        assert(nnz(I) == lp,'Second argument must be basis matrix.')

        cfp = p.coeffs(:,find(I));
    else
        cfp = [];
    end

    if all(~tf)
        % zero projection
        q = casadi.SX.zeros(lZ,1);
        e = p;
        return
    end
    
    % else:
    nT = z.nterm;

    % select coefficients
    Q = casadi.SX.zeros(nT,lp);
    Q(ii(tf),:) = cfp(idx(tf),:);

    % construct sparsity
    [ii,jj] = get_triplet(sparsity(z.coeffs));
    S = casadi.Sparsity.triplet(nT,lp,ii,floor(jj/lZ));

    % project onto template
    q = Q(find(S)); %#ok<FNDSB> 
end

% build projection error
if nargout < 3
    % no error
    e = casos.PS(0);
    return
end

% else
e = p - q'*z;

% e = casos.PS;
% 
% % make degree matrix unique
% [coeffs,degmat] = uniqueDeg(p.coeffs(idx(~tf),:), dgp(~tf,:));
% 
% % remove zero terms
% [coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);
% 
% % new polynomial
% e.coeffs = coeffs;
% e.degmat = degmat;
% e.indets = indets;
% e.matdim = size(p);

end
