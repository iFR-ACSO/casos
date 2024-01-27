function [coeffs,degmat] = uniqueDeg(coeffs,degmat,dim)
% Make degree matrix unique and return corresponding coefficients.
%
% This function ensures that the monomials are in graded REVERSE
% lexicographic order.

% sort by ascending degree
degsum = sum(degmat,2);

% make degree matrix unique
[degmat2,id,ic] = unique([degsum fliplr(degmat)],'rows','sorted');

if nargin < 3 || isempty(dim) || isequal('dim','all')
    % default or matrix operation over all elements
    % each row of the coefficient matrix corresponds to a term
    nt = size(coeffs,1);

    % sum coefficients for repeated terms
    summat = sparse(ic,1:nt,1,length(id),nt);
    coeffs = (summat*coeffs);

elseif dim == 1
    % each column of the coefficient matrix corresponds to a 
    % column of a coefficient per term, sorted by terms
    nc = size(coeffs,2);
    nt = size(degmat,1);

    % sum coefficient blocks for repeated terms
    summat = sparse(1:nc,kron(ic(:),ones(nc/nt,1)),1,nc,length(id)*nc/nt);
    coeffs = (coeffs*summat); %(summat*coeffs')';

elseif dim == 2
    % each row of the coefficient matrix corresponds to a 
    % row of a coefficient per term, sorted by row index
    nr = size(coeffs,1);
    nt = size(degmat,1);

    % sum coefficient rows for repeated terms
    summat = sparse(kron(ones(nr/nt,1),ic(:)),1:nr,1,length(id)*nr/nt,nr);
    coeffs = (summat*coeffs);

else
    error('Invalid dimension input.')
end

% reverse order of degrees
degmat = fliplr(degmat2(:,2:end));

end
