function z = scalar(vars,deg)
% Create a scalar monomial sparsity pattern.

if nargin == 0 || (nargin > 1 && isequal(deg, 0))
    % return zero-degree scalar
    z = casos.package.core.PolynomialSparsity.pattern(casadi.Sparsity.scalar);
    return
end

assert(is_indet(vars), 'First input must be indeterminate variables.')

p = casos.Indeterminates(vars);

if nargin < 2
    % return linear sparsity pattern
    degmat = speye(p.nvars);

else
    % return vector of monomials of desired degree
    assert(isvector(deg) && all(deg >= 0 & floor(deg) == ceil(deg)), 'Second input must be vector of nonnegative integers.')

    % enumerate monomials up to max(deg)
    r = nchoosek(p.nvars+max(deg),p.nvars); % total number of monomials
    M = nan(r,p.nvars);
    % iterate over ascending degrees
    l = 1;
    for i = unique(deg)
        [M,l] = degreemat(M,l,p.nvars,i);
    end

    % select degrees
    I = ismember(sum(M,2),deg);

    degmat = M(I,:);
end

% set output
z = build_monomials(degmat,p);

end

function [M,l] = degreemat(M, l0, m, d)
% Build degree matrix for m variables and degree d.
% Monomials will be in graded REVERSE lexicographic order.

    switch m
        case 1, M(l0,1) = d; l = l0+1;
        otherwise
            for j = d:-1:0
                [M,l] = degreemat(M,l0,m-1,j);
                M(l0:l-1,m) = d - j;
                l0 = l;
            end
    end
end
