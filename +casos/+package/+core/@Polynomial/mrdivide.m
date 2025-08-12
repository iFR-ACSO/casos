function c = mrdivide(p,b)
% Matrix right division with constant or symbolic matrix.
%
% Note: c = p/B is equivalent to solving the linear equation B*c = p.

assert(is_zerodegree(b),'Only division by constant or symbolic matrix possible.')

if isscalar(b)
    % fall back to scalar division
    c = rdivide(p,b);
    return
end

assert(~is_operator(b), 'Not allowed for operators.')

% input dimensions
szp = size(p);
szb = size(b);

if szb(2) ~= szp(2)
    % dimensions are compatible if size(b,2) == size(p,2)
    throw(casos.package.core.IncompatibleSizesError.matrix(p,b));
end

% else
c = b.new_poly;

% reshape coefficients
cfp = prepareMatrixOp(p.get_sparsity,p.coeffs,2);

% matrix division of coefficients
coeffs = cfp / p.new_coeff(b);

% output dimension (see MATLAB mrdivide)
sz = [szp(1) szb(1)];

% update sparsity pattern
[S,c.coeffs] = coeff_update(p.get_sparsity,coeffs,sz,2);

% set sparsity
c = set_sparsity(c,S);

end
