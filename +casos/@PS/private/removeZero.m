function [coeffs,degmat,indets] = removeZero(coeffs,degmat,indets)
% Remove degrees with zero coefficient, update degree matrix, 
% and remove unused variable.

nel = size(coeffs,2);

% remove degrees with zero coefficient
[coeffs,degmat] = removeCoeffs(coeffs,degmat);

% remove unused variables
[degmat, indets] = removeDegVar(degmat, indets);

% handle empty coeffient and/or degree matrix
if length(coeffs) < 1
    % no nonzero coefficients left
    coeffs = casadi.SX(1,nel);
    degmat = sparse(1,0);

elseif isempty(degmat)
    % constant polynomial
    coeffs = sum(coeffs,1);
    degmat = sparse(1,0);
end

end

function [coeffs,degmat] = removeCoeffs(coeffs,degmat)
% Remove terms with zero coefficient and update degree matrix.
    
    % nonzero coefficients with linear indices
    idx = find(sparsity(coeffs));
    [in,jn] = ind2sub(size(coeffs),idx);
    
    % identify all-zero rows
    [nr,~,ir] = unique(in);

    % length corresponds to number of nonzero rows
    I = 1:length(nr);
    
    % sparsity pattern without all-zero rows
    S = sparsity(casadi.SX(sparse(I(ir),jn,1)));
    % assign nonzero coefficients to sparsity pattern
    coeffs = casadi.SX(S,coeffs(idx));

    % remove corresponding degree matrix entries
    degmat = degmat(nr,:);
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
