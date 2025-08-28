function [S,coeffs] = coeff_int(obj,coeffs,x,range)
% Compute coefficient matrix for polynomial integral operations.

x = casos.Indeterminates(x);

% combine variables
[indets,dgf] = combineVar(obj.indets,x,obj.degmat);

% find position of integration variables
[tf,~] = ismember(indets,x);

% increase degrees
dgf(:,tf) = dgf(:,tf) + 1;

% compute scaling term for each monomial
scale = prod(dgf(:,tf),2);

if isempty(range)
    % indefinite integral
    degmat = dgf;

    if ~isa(coeffs,'casadi.Sparsity')
        % scale coefficients
        coeffs = scale.\coeffs;
    end

else
    % definite integral
    degmat = dgf(:,~tf);
    % remaining variables
    indets = indets(~tf);

    if ~isa(coeffs,'casadi.Sparsity')
        % integral(c_A*x^A, r1, r2) = prod_i c_A*[r2(i)^A(i) - r1(i)^A(i)]
        deg0 = dgf(:,tf)';
    
        % select range expressions
        [~,I] = sort(x);
        r = range(I,:);
    
        % compute exponents
        % prod_i [r2(i)^A(i) - r1(i)^A(i)]
        R = casadi_prod(r(:,2).^deg0 - r(:,1).^deg0,1);
    
        % coefficients of sum_A prod_i c_A*[r2(i)^A(i) - r1(i)^A(i)]
        coeffs = scale.\coeffs.*R';
    end

    % make degree matrix unique
    [coeffs,degmat] = uniqueDeg(coeffs, degmat);

    % remove zero terms
    [coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);
end

% return sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,size(obj));

end
