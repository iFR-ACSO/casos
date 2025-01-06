function [pts,P,V] = interpolation(S)
% Compute an interpolation basis based on the monomial sparsity pattern.
%
% The code is based on the approach(es) described in
%
% D. Papp and S. Yildiz. Sum-of-squares optimization without semidefinite 
% programming. Available at https://arxiv.org/abs/1712.01792.
%
% and
%
% A. Sommariva and M. Vianello, Computing approximate Fekete points by QR
% factorizations of Vandermonde matrices, Computers & Mathematics with
% Applications, 57 (2009), pp. 1324-1336. Available at
% https://doi.org/10.1016/j.camwa.2008.11.011.

if isempty(S)
    % empty sparsity pattern
    pts = casadi.DM(1,0);
    V = casadi.DM(0,0);
    P = casadi.DM(0,0);
    return

elseif is_zerodegree(S)
    % matrix sparsity pattern
    N = nnz(S);
    pts = casadi.DM(1,N);
    V = casadi.DM.eye(N);
    P = casadi.DM.eye(N);
    return
end

assert(isscalar(S), 'Not supported.')

% else
n = S.nvars;
d = S.maxdeg;

% aproximate Fekete points
% Chebyshev points of second type for each indeterminate variable
Pts = arrayfun(@(m) sin(pi*(-m:2:m)/(2*m)), d:(d+n-1), 'UniformOutput',false);
% generate possible combinations
pts = table2array(combinations(Pts{:}))';

% compute rectangular Vandermonde matrix
V = vandermonde(S,pts);

% compute vector of moments over domain [-1 1]
m = moments(S, [-1 1]);

% solve for quadrature weights
w = full(V')\full(m);

% detect nonzero weights
ind = find(w ~= 0);

assert(length(ind) == S.nterm, 'Something went wrong.')

% select full-rank submatrix of V
V = V(ind,:);

% keep only points with nonzero weights
pts = pts(:,ind);

% ortho-normalize basis
[P,~] = qr(V);

end

function m = moments(S,domain)
% Compute vector of moments.

    % integrate monomials over domain
    [~,m] = coeff_int(to_vector(S),casadi.DM.eye(S.nterm),S.indets,domain);
    
    % reshape to vector
    m = reshape(m,S.nterm,1);
end
