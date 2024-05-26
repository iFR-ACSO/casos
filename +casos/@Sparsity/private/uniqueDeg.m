function [coeffs,degmat,I,M] = uniqueDeg(coeffs,degmat)
% Make degree matrix unique and return corresponding coefficients.
%
% This function ensures that the monomials are in graded REVERSE
% lexicographic order.

nt = size(coeffs,1);
nv = size(coeffs,2);

if iscell(degmat)
    % undocumented: get indices
    [id,ic] = degmat{:};

else

% sort by ascending degree
degsum = sum(degmat,2);

% make degree matrix unique
[degmat2,id,ic] = unique([degsum fliplr(degmat)],'rows','sorted');

% reverse order of degrees
degmat = fliplr(degmat2(:,2:end));

end

if isa(coeffs,'casadi.Sparsity')
    % take union of repeated coefficients
    [ii,jj] = get_triplet(coeffs);
    coeffs = casadi.Sparsity.triplet(length(id),nv,ic(ii+1)-1,jj);

else
    % sum repeated coefficients
    summat = sparse(ic,1:nt,1,length(id),nt);
    coeffs = (summat*coeffs);
end

if nargout > 2
    % undocumented: return indices
    I = {id ic};
    M = sparse(ic,1:nt,1,length(id),nt);
end

end
