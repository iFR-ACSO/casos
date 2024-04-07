function z = monomials(vars,deg)
% Create monomial sparsity pattern.

assert(is_indet(vars), 'First input must be indeterminate variables.')

if deg == 0
    % Return constant one.
    degmat = sparse(1,0);
    indets = casos.Indeterminates;

else
    % Return vector of monomials of desired degree.
    assert(isvector(deg) && all(deg >= 0 & floor(deg) == ceil(deg)), 'Second input must be vector of nonnegative integers.')

    p = casos.Indeterminates(vars);

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
    indets = p;
end

% set output
z = build_monomials(degmat,indets);

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
