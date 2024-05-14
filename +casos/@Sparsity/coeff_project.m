function [S,coeffs] = coeff_project(obj,coeffs,S)
% Project polynomial coefficient matrix onto sparsity pattern.

assert(~isa(coeffs,'casadi.Sparsity'),'Notify the developers.')

% expand coefficients to match degrees
[S,S_coeffs,coeffs] = coeff_expand(S,obj,S.coeffs,coeffs);

if isa(coeffs,'casadi.DM') && ~is_regular(coeffs)
    % handle irregular coefficients (inf, nan)
    coeffs = extend_irregular(coeffs,+inf);
    coeffs = extend_irregular(coeffs,-inf);
    coeffs = extend_irregular(coeffs,nan);
end

% project onto coefficient sparsity pattern
coeffs = project(coeffs,S_coeffs);

% remove zero terms
[coeffs,S.degmat,S.indets] = removeZero(coeffs,S.degmat,S.indets);

% store coefficients
S = set_coefficients(S,coeffs);

end

function coeffs = extend_irregular(coeffs,num)
% If the coefficient matrix contains an irregular number, 
% extend that number over all terms.

    S_coeffs = sparsity(sparsify(coeffs == num));
    % detect columns (= entry of coefficients) that are irregular
    [~,j1] = get_triplet(S_coeffs);

    if isempty(j1), return; end % nothing to do

    % make entry irregular for all terms
    coeffs(:,unique(j1)+1) = num;
end