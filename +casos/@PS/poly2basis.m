function [q,z,e] = poly2basis(p,z,I)
% Project polynnomial p onto the span of monomials contained in z.

p = casos.PS(p);

if nargin < 2
    z = basis(p);
    I = true(size(p));
elseif islogical(z)
    I = z;
    z = basis(p,I);
elseif nargin < 3
    I = true(size(p));
end

% combine variables
[~,dgp,dgz] = combineVar(p.indets,z.indets,p.degmat,z.degmat);

% find degrees of z in p
[tf,ii] = ismember(dgp,dgz,'rows');
idx = 1:p.nterm;

if iscolumn(z) && all(I)
    % legacy code: project onto monomial vector
    assert(is_monom(z),'Second argument must be vector of monomials.')

    % select coefficients
    q = casadi.SX(numel(z),numel(p));

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
    lp = size(z,2);
    lZ = size(z,1);

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
        q = casadi.SX(lZ,1);
        e = p;
        return
    end
    
    % else:
    % select coefficients
    Q = casadi.SX(z.nterm,lp);
    Q(ii(tf),:) = cfp(idx(tf),:);

    % construct template
    S = z'*ones(lZ,1);

    % project onto template
    C = nonzeros(project(Q,sparsity(S.coeffs)));
    q = vertcat(C{:});
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
