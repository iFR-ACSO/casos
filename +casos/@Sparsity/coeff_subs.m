function [S,coeffs] = coeff_subs(S1,coeff1,x,S2,coeff2)
% Substitute indeterminates in coefficient matrix.

x = casos.Indeterminates(x);

% find position of indeterminates to substitute
[tf,xloc] = ismember(S1.indets,x);

if ~any(tf)
    % nothing to substitute
    S = S1;
    coeffs = coeff1;
    return
end

% else
S = casos.Sparsity;

% substitute b = sum_b c_b*y^b for x into a = sum_a*c_a*(x,y)^a
% yields c = sum_a c_a*b^a1*y^a2 with a = (a1,a2)
deg0 = S1.degmat(:, tf);

% select remaining coefficients and variables
indets = S1.indets(~tf);
degA = S1.degmat(:,~tf);

% select expressions to substitute with
[S2,coeff2] = coeff_subsref(S2,coeff2,xloc(tf),[1 nnz(tf)]);
% repeat to match number of lhs terms
[S2,coeff2] = coeff_repmat(S2,coeff2,S1.nterm,1);

% compute exponents b^a1
[Sb,coeffb] = coeff_power(S2,coeff2,deg0);
[Sb,coeffb] = coeff_prod(Sb,coeffb,2);

if is_zerodegree(S2)
    % substitution with scalar expression
    coeffb = repmat(coeffb,S1.numel,1);
    % coefficients of c = sum_a c_a*b^a1*y^a2
    coeffs = times_coeffs(coeff1, T(coeffb));
    % degree matrix of c 
    degmat = degA;

else
    % substitution with polynomial
    % combine variables
    [indets,degA,degB] = combineVar(indets,Sb.indets,degA,Sb.degmat);

    % compute coefficients for c
    % where b^a1 = sum_B c_B*y^B
    % hence c = sum_a sum_B c_a*c_B[a1] y^(a2+B) with a = (a1,a2)
    ntA = S1.nterm;
    ntB = Sb.nterm;

    % repeat matrices for double sum
    cfA = kron(coeff1,ones(ntB,1)); dgA = kron(degA,ones(ntB,1));
    cfB = kron(ones(ntA,1),coeffb); dgB = kron(ones(ntA,1),degB);

    % select c_B[a1]
    ii = 1:(ntA*ntB);
    jj = kron(1:ntA,ones(1,ntB));
    I = sub2ind(size(cfB),ii,jj);

    % combine coefficients and degrees
    if isa(cfB,'casadi.Sparsity')
        % combine sparsity patterns
        C = repmat(sub(cfB,I-1,casadi.Sparsity.dense(ntA*ntB,1)),1,S1.numel);
        coeffs = cfA * C;
    else
        % combine coefficient matrices
        C = reshape(cfB(I),ntA*ntB,1);
        coeffs = cfA .* C;
    end
    degmat = dgA + dgB;
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity
S.degmat = degmat;
S.indets = indets;
S.matdim = S1.matdim;
% store coefficients
S = set_coefficients(S,coeffs);

end
