function [coeffs,degmat] = uniqueDeg(coeffs,degmat,setOrder)
% Make degree matrix unique and return corresponding coefficients.

if nargin < 3
    setOrder = 'sorted';
end

nt = size(coeffs,1);

% make degree matrix unique
[degmat,id,ic] = unique(degmat,'rows',setOrder);

% sum repeated coefficients
summat = sparse(ic,1:nt,1,length(id),nt);
coeffs = sparsify(summat*coeffs);

% reverse order of degrees if sorted by unique
if isequal(setOrder,'sorted')
    coeffs = flipud(coeffs);
    degmat = flipud(degmat);
end

end
