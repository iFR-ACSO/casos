function Z = monomials(vars,deg)
% Create vector of monomials.

p = casos.PS(vars);

indets = p.indets;

if nargin == 1
    % Return vector of monomials in polynomial.

    % degrees are in graded reverse lexicographic order (canonic)
    degmat = p.degmat;

elseif deg == 0
    % Return constant one.
    degmat = sparse(1,0);
    indets = {};

else
    % Return vector of monomials of desired degree.
    assert(is_indet(p), 'First input must be vector of indeterminates.')
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
Z = build_monomials(degmat,indets);

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
