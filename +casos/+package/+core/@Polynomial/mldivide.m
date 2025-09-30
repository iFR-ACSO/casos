function c = mldivide(a,p)
% Matrix left division with constant or symbolic matrix.
%
% Note: c = A\p is equivalent to solving the linear equation A*c = p.

assert(is_zerodegree(a),'Only division by constant or symbolic matrix possible.')

if isscalar(a)
    % fall back to scalar division
    c = ldivide(a,p);
    return
end

assert(~is_operator(a), 'Not allowed for operators.')

% input dimensions
sza = size(a);
szp = size(p);

if szp(1) ~= sza(1)
    % dimensions are compatible if size(a,1) == size(p,1)
    throw(casos.package.core.IncompatibleSizesError.mdivide(a,p));
end

% else
c = p.new_poly;

% reshape coefficients
cfp = prepareMatrixOp(p.get_sparsity,p.coeffs,1);

% matrix division of coefficients
coeffs = p.new_coeff(a) \ cfp;

% output dimension (see MATLAB mldivide)
sz = [sza(2) szp(2)];

% update sparsity pattern
[S,c.coeffs] = coeff_update(p.get_sparsity,coeffs,sz,1);

% set sparsity
c = set_sparsity(c,S);

end
