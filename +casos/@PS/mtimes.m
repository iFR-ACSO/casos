function c = mtimes(a,b)
% Matrix multiplication of two polynomials.

a = casos.PS(a);
b = casos.PS(b);

if isscalar(a) || isscalar(b)
    % fall back to scalar multiplication
    c = times(a,b);
    return
end

% input dimensions
sza = size(a);
szb = size(b);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% dimensions are compatible if inner dimensions agree
assert(sza(2) == szb(1), errsz, sza, szb)

% TODO: handle or escape for other simple cases, e.g., scalar, constant
% matrix, single term etc.?

% else
c = casos.PS;

nta = a.nterm;
ntb = b.nterm;

% combine variables
[indets,dga,dgb] = combineVar(a.indets,b.indets,a.degmat,b.degmat);

% Vectorized code to compute coef matrix
% from multipoly
idx1 = reshape(1:sza(1)*sza(2),[sza(1),sza(2)]);
idx1 = repmat(idx1,[nta 1]);
idx1 = idx1(:);

idx2 = repmat(1:nta,[sza(1) 1]);
idx2 = repmat(idx2(:),[sza(2) 1]);

idx = sub2ind([sza(1)*sza(2) nta],idx1,idx2);
acoef = a.coeffs;
acoefcol = full(acoef');
acoefcol = reshape(acoefcol(idx),[nta*sza(1) sza(2)]);

bcoef = b.coeffs;
bcoefcol = reshape(bcoef',[szb(1) szb(2)*ntb]);

tempcoef = acoefcol*bcoefcol;

idx1 = reshape(1:sza(1)*nta,[sza(1),nta]);
idx1 = repmat(idx1,[ntb*szb(2) 1]);
idx1 = idx1(:);

idx2 = repmat(1:(ntb*szb(2)),[sza(1) 1]);
idx2 = repmat(idx2(:),[nta 1]);

idx = sub2ind([nta*sza(1) ntb*szb(2)],idx1,idx2);
coeffs = reshape(tempcoef(idx),sza(1)*szb(2),nta*ntb)';

% (sum_a c_a*x^a)*(sum_b c_b*x^b) = (sum_a sum_b (c_a*c_b)*(x^a*x^b)
degmat = kron(dga,ones(ntb,1)) +  kron(ones(nta,1),dgb);

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new polynomial
c.coeffs = coeffs;
c.degmat = degmat;
c.indets = indets;
c.matdim = [sza(1) szb(2)];

end