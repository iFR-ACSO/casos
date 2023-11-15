function [coeffs,degmat] = uniqueDeg(coeffs,degmat)
% Make degree matrix unique and return corresponding coefficients.
%
% This function ensures that the monomials are in graded REVERSE
% lexicographic order.

nt = size(coeffs,1);

% make degree matrix unique
[degmat2,id,ic] = unique(fliplr(degmat),'rows','sorted');

% sum repeated coefficients
summat = sparse(ic,1:nt,1,length(id),nt);
coeffs = sparsify(summat*coeffs);

% reverse order of degrees
degmat = fliplr(degmat2);

end
