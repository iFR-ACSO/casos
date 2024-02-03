function c = ptaylor(a,x,b,deg)
% Polynomial Taylor expansion.

assert(nargin == 4, 'Incorrect number of inputs (got: %d, expected 4).',nargin)

a = casos.PS(a);
b = casos.PS(b);

assert(is_indet(x),'Second argument must be vector of indeterminate variables.')

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
c = casos.PS;

% shift to operating point
% TODO: use internal operation
A = subs(a,x,x-b);

% get degrees
D = sum(A.degmat,2);

% terms with exceeding degrees
coeffs = subsasgn(A.coeffs, substruct('()',{find(D > deg) ':'}), 0);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,A.degmat,A.indets);

% new polynomial
c.coeffs = coeffs;
c.degmat = degmat;
c.indets = indets;
c.matdim = A.matdim;

end
