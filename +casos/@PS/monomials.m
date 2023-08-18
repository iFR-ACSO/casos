function Z = monomials(vars,deg)
% Create vector of monomials.

p = casos.PS(vars);

indets = p.indets;

if nargin == 1
    % Return vector of monomials in polynomial.

    % TODO: can we assume that degrees are already sorted?
    % from multipoly:
    % Flipping/sorting to return output in lexicographic order
    % (sorted by degree then by alphabetical order)
    degmatsort = sortrows([sum(p.degmat,2) fliplr(p.degmat)]);
    degmat = fliplr(degmatsort(:,2:end));

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
    % iterate over degrees
    l = 1;
    for i = deg
        [M,l] = degreemat(M,l,p.nvars,i);
    end

    % select degrees
    I = ismember(sum(M,2),deg);

    degmat = M(I,:);
end

% number of monomials
nt = size(degmat,1);

% set output
Z = casos.PS;
Z.coeffs = casadi.SX.eye(nt);
Z.degmat = sparse(degmat);
Z.indets = indets;
Z.matdim = [nt 1];

end

function [M,l] = degreemat(M, l0, m, d)
% Build degree matrix for m variables and degree d.

    switch m
        case 1, M(l0,end) = d; l = l0+1;
        otherwise
            for j = 0:d
                [M,l] = degreemat(M,l0,m-1,j);
                M(l0:l-1,end-m+1) = d - j;
                l0 = l;
            end
    end
end
