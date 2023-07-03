function c = plus(a,b)
% Add two polynomials.

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

% dimensions are compatible if equal or one summand is row/column
assert(all(I | I1), errsz, sza, szb)

% dimensions of sum
sz = max(sza,szb);

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % addition with empty polynomial is empty
    sz(I0) = 0;

    c = casos.PS.zeros(sz);
    return
end

% else
c = casos.PS;

nta = a.nterm; 
ntb = b.nterm;

% reshape to output dimensions
cfa = reshape(repmat(a.coeffs,sz./sza),nta,prod(sz));
cfb = reshape(repmat(b.coeffs,sz./szb),ntb,prod(sz));

% combine variables
[indets,dga,dgb] = combineVar(a.indets,b.indets,a.degmat,b.degmat);

% make degree matrix unique
[coeffs,degmat] = uniqueDeg([cfa;cfb], [dga;dgb]);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new polynomial
c.coeffs = flipud(coeffs);
c.degmat = flipud(degmat);
c.indets = indets;
c.matdim = sz;

end