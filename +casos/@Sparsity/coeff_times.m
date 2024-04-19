function [S,coeffs] = coeff_times(S1,S2,cfa,cfb)
% Multiply (intersect) two polynomial coefficient matrices.

S = casos.Sparsity;

nta = S1.nterm;
ntb = S2.nterm;

% combine variables
[indets,dga,dgb] = combineVar(S1.indets,S2.indets,S1.degmat,S2.degmat);

% (sum_a c_a*x^a).*(sum_b c_b*x^b) = sum_a sum_b (c_a.*c_b)*(x^a*x^b)
coeffs = times_coeffs(kron(cfa,ones(ntb,1)), kron(ones(nta,1),cfb));
degmat = times_degmat(kron(dga,ones(ntb,1)), kron(ones(nta,1),dgb));

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S.degmat = degmat;
S.indets = indets;
S.matdim = size(S1);
% store coefficients
S = set_coefficients(S,coeffs);

end

function C = times_coeffs(A,B)
    % Multiply coefficient matrices.

    if isa(A,'casadi.Sparsity')
        % intersect sparsity
        C = intersect(A,B);

    else
        % multiply coefficients element-wise
        C = times(A,B);
    end
end

function D = times_degmat(A,B)
    % Multiply degree matrices
    D = plus(A,B);
end
