function c = kron(a,b)
% Compute Kronecker product of two polynomial matrices.

a = casos.PS(a);
b = casos.PS(b);

% input dimensions
sza = size(a);
szb = size(b);

% dimensions of product
sz = sza.*szb;

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % product with empty polynomial is empty

    c = casos.PS.zeros(sz);
    return
end

% else
c = casos.PS;

% coefficients
cfa = a.coeffs;
cfb = b.coeffs;

% number of terms & elements
[nta,nea] = size(cfa);
[ntb,neb] = size(cfb);

% combine variables
[indets,~,ic] = unique([a.indets b.indets]);

nv = length(indets);

% get nonzero degrees
[ia,ja,da] = find(a.degmat);
[ib,jb,db] = find(b.degmat);

% extend degree matrices to combined variables
dga = sparse(ia,ic(ja),da,nta,nv);
dgb = sparse(ib,ic(a.nvars+jb),db,ntb,nv);

% reshape coefficients
[Ia,Ib] = ndgrid(1:nta,1:ntb);
ja = reshape(kron(reshape(1:nea,sza),ones(szb)),[],1);
jb = reshape(kron(ones(sza),reshape(1:neb,szb)),[],1);

degmat = dga(Ia(:),:) + dgb(Ib(:),:);
coeffs = cfa(Ia(:),ja(:)).*cfb(Ib(:),jb(:));

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs,degmat);

% new polynomial
c.coeffs = coeffs;
c.degmat = degmat;
c.indets = indets;
c.matdim = sz;

end