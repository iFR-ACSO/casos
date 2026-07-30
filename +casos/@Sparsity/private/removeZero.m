% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [coeffs,degmat,indets] = removeZero(coeffs,degmat,indets)
% Remove degrees with zero coefficient, update degree matrix, 
% and remove unused variable.
%
% Arguments:
%
%   coeffs                         matrix of which the rows are vec(c_a)
%   degmat                         matrix of which the rows are [a1 ... aN]
%   indets                         variables {x1,...,xN} 
%
% The number of rows in coeffs and degmat equals the number of terms; after
% evaluation, this equals the number of terms with nonzero coefficients.
% The number of columns in coeffs equals the number of elements and remains
% unchanged.
% The number of columns in degmat equals the number of indeterminate
% variables (N); if the indets argument is given, indeterminates with
% all-zero degrees are removed; otherwise, this remains unchanged.

% remove degrees with zero coefficient
[coeffs,degmat] = removeCoeffs(coeffs,degmat);

if nargin > 2
% remove unused variables
[degmat, indets] = removeDegVar(degmat, indets);
end

if isempty(degmat) || size(coeffs,1) < 1
    % constant polynomial
    coeffs = sum1(coeffs);
end

end

function [coeffs,degmat] = removeCoeffs(coeffs,degmat)
% Remove terms with zero coefficient and update degree matrix.
    
    if isa(coeffs,'casadi.Sparsity')
        % nonzero row and column indices (triplet)
        [ii,jj] = ind2sub(size(coeffs),find(coeffs));

        % sparsity pattern without all-sparse rows
        [coeffs,nr] = nonzeroPattern(size(coeffs),ii,jj);
        
    else
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
    
        % sparsity pattern without all-sparse rows
        [S,nr] = nonzeroPattern(size(coeffs),ii,jj);

        % assign nonzero coefficients to sparsity pattern
        coeffs = project(coeffs(nr,:),S);
    end

    if isempty(ii)
        % all-sparse coefficients
        % return row vectors for single (constant) monomial term
        degmat = sparse(1,size(degmat,2));

    else
        % remove corresponding degree matrix entries
        degmat = degmat(nr,:);
    end
end

function [S,nr] = nonzeroPattern(sz,ii,jj)
% Return sparsity pattern without all-sparse rows

    % identify all-zero rows (zero coefficient)
    [nr,~,ir] = unique(ii);

    % length corresponds to number of nonzero rows
    I = 1:length(nr);
    
    % sparsity pattern without all-sparse rows (note: 0-based index)
    S = casadi.Sparsity.triplet(length(nr),sz(2),I(ir)-1,jj-1);
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
