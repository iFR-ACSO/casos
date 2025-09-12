function [coeffs,Si,So] = removeZero(coeffs,Si,So)
% Remove zero coefficients.

if isa(coeffs,'casadi.Sparsity')
    % nonzero row and column indices (triplet)
    [ii,jj] = get_triplet(coeffs);  % CasADi has 0-index

    % sparsity pattern without all-sparse rows or columns
    [coeffs,I,J] = nonzeroPattern(ii,jj);

else
    % nonzero row and column indices (triplet)
    [ii,jj] = get_triplet(coeffs);  % CasADi has 0-index

    % sparsity pattern without all-sparse rows or columns
    [S,I,J] = nonzeroPattern(ii,jj);

    % assign nonzero coefficients
    coeffs = project(coeffs(I,J),S);
end

% remove corresponding monomial entries
Si = coeff_setnz(Si,J);
So = coeff_setnz(So,I);

end

function [S,nr,nc] = nonzeroPattern(ii,jj)
% Return sparsity pattern without all-sparse rows or columns.

    % identify all-zero rows
    [nr,~,ir] = unique(ii+1);   % return 1-index
    % identify all-zero columns
    [nc,~,ic] = unique(jj+1);   % return 1-index

    % size corresponds to number of nonzero rows and columns
    I = 1:length(nr);
    J = 1:length(nc);

    % sparsity pattern with all-sparse rows (0-based index)
    S = casadi.Sparsity.triplet(length(nr),length(nc),I(ir)-1,J(ic)-1);
end
