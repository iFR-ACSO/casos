function p = cleanpoly(p,tol,deg)
% Remove coefficients from nonsymbolic polynomials.

% check inputs
if nargin<3
    deg = [];
end

if ~isempty(tol)
    % remove coefficients below tolerance
    assert(isdouble(tol) && isscalar(tol) && tol >= 0, 'Tolerance must be a positive scalar or is empty.')

    p = remove_coeffs(p,tol);
end

if ~isempty(deg)
    % retain terms with certain degrees
    assert(isa(deg,'double') &&  all(floor(deg)==ceil(deg)) && all(deg>=0), 'Third input must be a vector of non-negative integers.')

    % restrict polynomial sparsity to given degrees
    S = restrict_terms(sparsity(p), deg);

    % project onto restricted sparsity pattern
    p = project(p,S);
end

end
