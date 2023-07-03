function [coeffs,degmat] = uniqueDeg(coeffs,degmat)
% Make degree matrix unique and return corresponding coefficients.

nt = size(coeffs,1);

% make degree matrix unique
[degmat,id,ic] = unique(degmat,'rows','sorted');

% sum repeated coefficients
summat = sparse(ic,1:nt,1,length(id),nt);
coeffs = sparsify(summat*coeffs);

end
