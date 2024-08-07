function [S,coeffs] = coeff_transpose(obj,coeffs)
% Transpose of a polynomial coefficient matrix.

S = casos.Sparsity;

% dimensions
sza = size(obj);
szb = [sza(2) sza(1)];

if ~isvector(obj)
    % transpose matrix

    % indices
    [ib,jb] = ind2sub(szb,1:prod(szb));
    
    % transpose
    idx = sub2ind(sza,jb,ib);
    
    if isa(coeffs,'casadi.Sparsity')
        % subindex into sparsity pattern
        % NOTE: Casadi has zero-based index
        coeffs = sub(coeffs,0:(obj.nterm-1),idx-1);
    
    else
        % subscript into coefficient matrix
        coeffs = coeffs(:,idx);
    end
end
% nothing to do for a vector

% new sparsity pattern
S.degmat = obj.degmat;
S.indets = obj.indets;
S.matdim = szb;
% store coefficients
S = set_coefficients(S,coeffs);

end
