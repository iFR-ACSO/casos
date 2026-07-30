% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [c,S] = poly2basis(obj,S)
% Return a vector of nonzero coordinates for a given basis (sparsity).
% If no basis is given, the polynomial's sparsity pattern is used.

if nargin < 2
    % return nonzero coefficients (below)
    S = sparsity(obj);
    coeffs = obj.coeffs;

elseif ~isscalar(obj) && ~islogical(S) ...
         && ~(all(size(obj) == 0) && isempty(S)) && ~check_sz_equal(obj,S)
    % dimension mismatch
    throw(casos.package.core.IncompatibleSizesError.other(obj,S));

elseif isempty(obj) || (isscalar(obj) && isempty(S)) || (islogical(S) && all(~S))
    % empty polynomial, selection, or projection
    assert(~isempty(obj) || isempty(S),'Cannot project empty polynomial.')

    c = obj.new_coeff(0,1);
    S = casos.Sparsity(0,0);
    return

elseif islogical(S)
    % nonzeros and basis of subreference
    [i,j] = coeff_triplet(obj.get_sparsity);
    tf = S(j+1);
    idx = sub2ind(size(obj.coeffs),i(tf)+1,j(tf)+1);
    c = reshape(obj.coeffs(idx),nnz(tf),1);
    S = basis(obj,S);
    return

elseif isscalar(obj) && ~isscalar(S)
    % repeat scalar inputs before projection
    [sp_rep,coeffs] = coeff_repmat(obj.get_sparsity,obj.coeffs,size(S));
    [S,coeffs] = coeff_project(sp_rep,coeffs,S,true);

else
    % project onto given sparsity pattern.
    [S,coeffs] = coeff_project(obj.sparsity,obj.coeffs,S,true);
end

% select nonzero coefficients
c = reshape(coeffs(coeff_find(S)),nnz(S),1);

end
