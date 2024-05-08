function c = mtimes(a,b)
% Matrix multiplication of two polynomials.

if isscalar(a) || isscalar(b)
    % fall back to scalar multiplication
    c = times(a,b);
    return

elseif ~check_sz_mtimes(a,b)
    % dimensions are compatible if size(a,2) == size(b,1)
    throw(casos.package.core.IncompatibleSizesError.matrix(a,b));
end

% else
c = a.new_poly;

nta = a.nterm;
ntb = b.nterm;

sza = size(a);
szb = size(b);

if isrow(a) && isvector(b)
    % inner vector product
    % (sum_a c_a'*x^a)*(sum_b c_b*x^b) = sum_a sum_b (c_a'*c_b)*(x^a*x^b)
    coeffs = reshape(a.coeffs*b.coeffs',nta*ntb,1);

else
% Vectorized code to compute coef matrix
% from multipoly
idx1 = reshape(1:sza(1)*sza(2),[sza(1),sza(2)]);
idx1 = repmat(idx1,[nta 1]);
idx1 = idx1(:);

idx2 = repmat(1:nta,[sza(1) 1]);
idx2 = repmat(idx2(:),[sza(2) 1]);

idx = sub2ind([sza(1)*sza(2) nta],idx1,idx2);
acoef = a.coeffs;
acoefcol = (acoef');
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
end

% update coefficients and degrees
[S,c.coeffs] = coeff_mtimes(a.get_sparsity,b.get_sparsity,coeffs);

% set sparsity
c = set_sparsity(c,S);

end
