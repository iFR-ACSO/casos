function [c,pts] = poly2interpol(obj,pts)
% Return a vector of coordinates for a given interpolation.
% If no interpolation is given, the interpolation basis is computed based
% on the polynomial's sparsity pattern.

assert(isscalar(obj), 'Not supported.')

if nargin < 2
    % compute interpolation basis
    pts = interpolation(obj.poly_sparsity);
end

nt = obj.nterm;
nv = obj.nvars;

% convert polynomial to function
X = casadi.SX.sym('x',nv,1);
% substitute indeterminate variables
[~,px] = coeff_subs(...
            obj.poly_sparsity, ...
            obj.coeffs, ...
            obj.indeterminates, ...
            casos.Sparsity.dense(nv,1), ...
            X' ...
);
% to function
pfun = casadi.Function('p',{X},{sum(px,2)},struct('allow_free',true));
% pfun = to_function(obj,struct('allow_free',true));

% evaluate on interpolation
c = pfun(obj.new_coeff(pts))';

end
