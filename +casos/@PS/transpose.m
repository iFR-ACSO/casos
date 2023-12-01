function b = transpose(a)
% Non-conjugate transpose of a polynomial matrix.

b = casos.PS;

% dimensions
sza = size(a);
szb = [sza(2) sza(1)];

% indices
[ib,jb] = ind2sub(szb,1:prod(szb));

% transpose
idx = sub2ind(sza,jb,ib);

b.coeffs = a.coeffs(:,idx);
b.degmat = a.degmat;
b.indets = a.indets;
b.matdim = szb;

end