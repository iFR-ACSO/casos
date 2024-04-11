function [coeffs,degmat,indets] = removeZero(coeffs,degmat,indets)
% Remove degrees with zero coefficient, update degree matrix, 
% and remove unused variable.

nel = size(coeffs,2);

if ~isa(coeffs,'casadi.Sparsity')
% remove degrees with zero coefficient
[coeffs,degmat] = removeCoeffs(coeffs,degmat);
end

if nargin > 2
% remove unused variables
[degmat, indets] = removeDegVar(degmat, indets);
end

% handle empty coeffient and/or degree matrix
if length(coeffs) < 1
    % no nonzero coefficients left
    error('Notify the developers')
    coeffs = casadi.SX(1,nel);
    degmat = sparse(1,0);

elseif isempty(degmat)
    % constant polynomial
    coeffs = sum1(coeffs);
    degmat = sparse(1,0);
end

end

function [coeffs,degmat] = removeCoeffs(coeffs,degmat)
% Remove terms with zero coefficient and update degree matrix.
    
    % sparsity pattern of coefficients
    S1 = sparsity(coeffs);
    % sparsity pattern without zeros
    S0 = sparsity(sparsify(coeffs));

    % nonzero coefficients with linear indices
    [i0,j0] = ind2sub(size(coeffs),find(S0));   % sparse zeros
    [~, j1] = ind2sub(size(coeffs),find(S1));   % full zeros
    
    % detect all-zero columns (full zero)
    iz = ~ismember(j1,j0);

    % assign full zeros to first not all-sparse row
    ii = [i0 repmat(max([min(i0) 1]),1,nnz(iz))];
    jj = [j0 j1(iz)];
    
    % identify all-zero rows (zero coefficient)
    [nr,~,ir] = unique(ii);

    % length corresponds to number of nonzero rows
    I = 1:length(nr);
    
    % sparsity pattern without all-sparse rows (note: 0-based index)
    S = casadi.Sparsity.triplet(length(nr),size(coeffs,2),I(ir)-1,jj-1);
    % assign nonzero coefficients to sparsity pattern
    coeffs = project(coeffs(nr,:),S);

    if isempty(i0)
        % only zero coefficients
        degmat = sparse(1,size(degmat,2));
    else
        % remove corresponding degree matrix entries
        degmat = degmat(nr,:);
    end
end

function [degmat,indets] = removeDegVar(degmat,indets)
% Remove indeterminates with zero degree.

    % row indices of nonzero degrees
    [~,jd] = find(degmat);

    % identify all-zero columns
    nc = unique(jd);
    % remove all-zero columns
    degmat = degmat(:,nc);
    
    % remove corresponding indeterminates
    indets = indets(nc);
end
