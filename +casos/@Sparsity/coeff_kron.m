function [S,coeffs] = coeff_kron(S1,S2,coeff1,coeff2)
% Kronecker product of polynomial coefficient matrices.

S = casos.Sparsity;

% input dimensions
sza = size(S1);
szb = size(S2);

% number of terms & elements
[nta,nea] = size(coeff1);
[ntb,neb] = size(coeff2);

% combine variables
[indets,dga,dgb] = combineVar(S1.indets,S2.indets,S1.degmat,S2.degmat);

% reshape coefficients and degree matrix
% TODO: restrict to nonzero coefficients
[Ia,Ib] = ndgrid(1:nta,1:ntb);
ja = reshape(kron(reshape(1:nea,sza),ones(szb)),[],1);
jb = reshape(kron(ones(sza),reshape(1:neb,szb)),[],1);

if isa(coeff1,'casadi.Sparsity')
    % intersect coefficient patterns
    coeffs = sub(coeff1,Ia(:)-1,ja(:)-1)*sub(coeff2,Ib(:)-1,jb(:)-1); % NOTE: Casadi has 0-based index
else
    % multiply coefficient matrices
    coeffs = coeff1(Ia(:),ja(:)).*coeff2(Ib(:),jb(:));
end
% combine degrees
degmat = dga(Ia(:),:) + dgb(Ib(:),:);

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs,degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S.degmat = degmat;
S.indets = indets;
S.matdim = sza.*szb;
% store coefficients
S = set_coefficients(S,coeffs);

end
