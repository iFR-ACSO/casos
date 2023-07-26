function b = transpose(a)
% Non-conjugate transpose of a polynomial matrix.

b = casos.PS;

% dimensions
sza = size(a);
szb = [sza(2) sza(1)];

% indices
[ia,ja] = ind2sub(sza,1:prod(sza));

% transpose
idx = sub2ind(szb,ja,ia);

b.coeffs = a.coeffs(:,idx);
b.degmat = a.degmat;
b.indets = a.indets;
b.matdim = szb;

end