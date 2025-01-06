function V = vandermonde(S,pts)
% Return Vandermonde matrix with given discretization points.
%
% The Vandermonde matrix can be used to solve the interpolation problem
%
%   V(pts)*c = f(pts)
%
% where c is a vector of (unknown) coefficients and f is the function to be
% interpolated.
%
% Specifically, the Vandermonde matrix is defined as
%
%       | p1^a1 ... p1^aN |
%   V = |   :         :   |
%       | pM^a1 ... pM^aN |
%
% where pts = (p1,...,pM) are the interpolation points and (a1,...,aN) are
% the degrees of the monomial basis in multi-index notation.

if isempty(S)
    % empty sparsity pattern
    V = casadi.DM(0,0);
    return

elseif is_zerodegree(S)
    % matrix sparsity pattern
    V = casadi.DM.eye(nnz(S));
    return
end

assert(isscalar(S), 'Not supported.')

nt = S.nterm;
nv = S.nvars;

% convert monomial basis to (casadi) function
X = casadi.SX.sym('x',nv,1);
% substitute indeterminate variables
[~,px] = coeff_subs(...
            to_vector(S), ...
            casadi.DM.eye(nt), ...
            S.indets, ...
            casos.Sparsity.dense(nv,1), ...
            X' ...
);
% to function
pfun = casadi.Function('p',{X},{reshape(px,nt,1)});

% evaluate monomials in given points
V = pfun(pts)';

end
