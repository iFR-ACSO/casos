function c = times(a,b)
% Element-wise multiplication of two polynomial matrices.

a = casos.PS(a);
b = casos.PS(b);

% input dimensions
sza = size(a);
szb = size(b);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% compare dimensions
I = (sza == szb);
% find zero dimension
I0 = (sza == 0) | (szb == 0);
% find one dimension
I1 = (sza == 1) | (szb == 1);

% dimensions are compatible if equal or one factor is row/column
assert(all(I | I1), errsz, size2str(sza), size2str(szb))

% dimensions of element-wise product
sz = max(sza,szb);

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % element-wise multiplication with empty polynomial is empty
    sz(I0) = 0;

    c = casos.PS.zeros(sz);
    return
end
% TODO: handle or escape for other simple cases, e.g., scalar, constant
% matrix, single term etc.?

% else
c = casos.PS;

nta = a.nterm; 
ntb = b.nterm;

% reshape to output dimensions
cfa = reshape(repmat(a.coeffs,sz./sza),nta,prod(sz));
cfb = reshape(repmat(b.coeffs,sz./szb),ntb,prod(sz));

% combine variables
[indets,dga,dgb] = combineVar(a.indets,b.indets,a.degmat,b.degmat);

% (sum_a c_a*x^a).*(sum_b c_b*x^b) = sum_a sum_b (c_a.*c_b)*(x^a*x^b)
coeffs = kron(cfa,ones(ntb,1)) .* kron(ones(nta,1),cfb);
degmat = kron(dga,ones(ntb,1)) +  kron(ones(nta,1),dgb);

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new polynomial
c.coeffs = coeffs;
c.degmat = degmat;
c.indets = indets;
c.matdim = sz;

end