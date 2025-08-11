function [tf,I] = is_monom(obj)
% Check if polynomial sparsity pattern is a vector of monomial terms.

[ia,ja] = get_triplet(obj.coeffs);

% pattern is a vector of monomial terms 
% if all coefficients have exactly one nonzero entry
tf = issorted(ja,'strictascend');

% return order of monomials in vector
I = ia + 1; % NOTE: Casadi has zero-based index

end
