function c = subs(a,x,b)
% Polynomial substitution of indeterminates x in a by expression b.

a = casos.PS(a);
b = casos.PS(b);

assert(is_indet(x),'Second argument must be vector of indeterminate variables.')
assert(length(x) == length(b),'Second and third argument have incompatible sizes.')

[tf,xloc] = ismember(a.indets,x.indets);

if ~any(tf)
    % nothing to substitute
    c = a;
    return
end

% else:
c = casos.PS;

% select expressions to substitute with
b = x.coeffs(xloc(tf),:)*b; %TODO: internal operation

% substitute b = sum_b c_b*y^b for x into a = sum_a*c_a*(x,y)^a
% yields c = sum_a c_a*b^a1*y^a2 with a = (a1,a2)
deg0 = a.degmat(:, tf);

% select remaining coefficients and variables
indets = a.indets(~tf);
degA = a.degmat(:,~tf);

if is_zerodegree(b) 
    % substitution with scalar expression
    b = repmat(reshape(b.coeffs,1,length(b)),a.nterm,1);
    
    % compute exponents b^a1
    B = prod(b.^deg0,2);

    % coefficients of c = sum_a c_a*b^a1*y^a2
    coeffs = a.coeffs.*B;
    % degree matrix of c 
    degmat = degA;

else
    % substitution with polynomial
    assert(length(x) == 1,'Not supported yet.')

    % compute exponents of b'
    B = (b').^deg0;

    % combine variables
    [indets,degA,degB] = combineVar(indets,B.indets,degA,B.degmat);

    % compute coefficients for c
    % where b^a1 = sum_B c_B*y^B
    % hence c = sum_a sum_B c_a*c_B[a1] y^(a2+B) with a = (a1,a2)
    error('Not implemented.')
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new polynomial
c.coeffs = coeffs;
c.degmat = degmat;
c.indets = indets;
c.matdim = a.matdim;

end

function B = prod(A,dim)
% Product of matrix elements.

    if nargin < 2 && isrow(A)
        dim = 2;
    elseif nargin < 2
        dim = 1;
    end

    % size of result
    sz = size(A);
    sz(dim) = 1;
    % prepare result
    B = ones(sz);
    % select elements to multiply
    L.type = '()'; L.subs = {':',':'};

    for i = 1:size(A,dim)
        % select elements
        L.subs{dim} = i;
        B = B.*subsref(A,L);
    end
end