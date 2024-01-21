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
    B = sx_prod(b.^deg0,2);

        % retain constant monomials (work around)
    B(find(all(deg0==0,2))) = 1; %#ok<FNDSB> 
    
    % coefficients of c = sum_a c_a*b^a1*y^a2
    coeffs = a.coeffs.*B;
    % degree matrix of c 
    degmat = degA;

else
    % substitution with polynomial
    b = reshape(b,1,length(b));

    % compute exponents of b'
    B = prod(b.^deg0,2);

    % combine variables
    [indets,degA,degB] = combineVar(indets,B.indets,degA,B.degmat);

    % compute coefficients for c
    % where b^a1 = sum_B c_B*y^B
    % hence c = sum_a sum_B c_a*c_B[a1] y^(a2+B) with a = (a1,a2)
    ntA = a.nterm;
    ntB = B.nterm;

    % repeat matrices for double sum
    cfA = kron(a.coeffs,ones(ntB,1)); dgA = kron(degA,ones(ntB,1));
    cfB = kron(ones(ntA,1),B.coeffs); dgB = kron(ones(ntA,1),degB);

    % select c_B[a1]
    idx = kron(1:ntA,ones(1,ntB));
    I = sub2ind(size(cfB),1:(ntA*ntB),idx);

    % combine coefficients and degrees
    coeffs = cfA .* cfB(I');
    degmat = dgA + dgB;
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
