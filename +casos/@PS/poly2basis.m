function [q,z,e] = poly2basis(p,z)
% Project polynnomial p onto the span of monomials contained in z.

p = casos.PS(p);

if nargin < 2
    z = monomials(p);
else
    assert(is_monom(z),'Second argument must be vector of monomials.')
end

% combine variables
[indets,dgp,dgz] = combineVar(p.indets,z.indets,p.degmat,z.degmat);

% find degrees of z in p
[tf,I] = ismember(dgp,dgz,'rows');
idx = 1:p.nterm;

% select coefficients
q = casadi.SX(numel(z),numel(p));
q(I(tf),:) = p.coeffs(idx(tf),:);
% reorder to match monomials
q = z.coeffs * q;

% build projection error
if all(tf)
    % no error
    e = casos.PS(0);
    return
end

% else
e = casos.PS;

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(p.coeffs(idx(~tf),:), dgp(~tf,:));

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new polynomial
e.coeffs = coeffs;
e.degmat = degmat;
e.indets = indets;
e.matdim = size(p);

end
