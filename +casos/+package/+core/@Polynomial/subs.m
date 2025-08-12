function c = subs(a,x,b)
% Polynomial substitution of indeterminates x in a by expression b.

assert(~is_operator(a),'First argument must not be an operator.')
assert(is_indet(x),'Second argument must be vector of indeterminate variables.')
assert(~is_operator(b),'Third argument must not be an operator.')
        
% check dimensions
if isscalar(b)
    % replace all indeterminates by same expression
    b = repmat(b,size(x));

elseif (size(x,1) == size(b,1) && iscolumn(x) && size(b,2) > 1)
    % repeated substitution -- not supported
    error('Repeated substitution not supported, use to_function(a) instead.')

else
    assert(numel(x) == numel(b),'Second and third argument have incompatible sizes.')
end

c = a.new_poly;

% substitute indeterminates in coefficient matrix
[S,c.coeffs] = coeff_substitute(a.get_sparsity,a.coeffs,x,b.get_sparsity,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
