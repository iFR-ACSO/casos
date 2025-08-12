function c = ptaylor(a,x,b,deg)
% Polynomial Taylor expansion.

assert(nargin == 4, 'Incorrect number of inputs (got: %d, expected 4).',nargin)
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

% else:
c = a.new_poly;

% shift to operating point
y = x + b;
% TODO: use internal operation
[S,coeffs] = coeff_subs(a.get_sparsity,a.coeffs,x,y.get_sparsity,y.coeffs);

% remove terms with exceeding degrees
[S,c.coeffs] = coeff_project(S,coeffs,restrict_terms(S,0:deg));

% set sparsity
c = set_sparsity(c,S);

end
