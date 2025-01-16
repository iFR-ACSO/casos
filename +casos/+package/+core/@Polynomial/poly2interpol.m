function [c,pts,S] = poly2interpol(obj,pts,S)
% Return a vector of coordinates for a given interpolation.
% If no interpolation is given, the interpolation basis is computed based
% on the polynomial's sparsity pattern.

if nargin < 3
    % use sparsity pattern
    S = obj.poly_sparsity;
end
if nargin < 2
    % compute interpolation basis
    pts = interpolation(S);
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
pfun = casadi.Function('p',{X},{px'},struct('allow_free',true));
% pfun = to_function(obj,struct('allow_free',true));

% evaluate on interpolation
c = pfun(obj.new_coeff(pts))';

if ~isscalar(obj)
    % multidimensional polynonial [p1; ...; pn]
    % return [p1(t1) ... p1(tL), ..., pn(t1) ... pn(tN)]
    % where (L, ..., N) corresponds to the number of terms in p1, ..., pn,
    % respectively.
    nts = get_nterm(S)';

    % iterate elements in matrix c
    I = 1:numel(c);
    % column indices
    [ii,jj] = ind2sub(size(c),I);
    % rows correspond to number of evaluation points
    tf = (ii <= nts(jj));
    % remove additional points
    I(~tf) = [];

    % return
    c = c(I)';
end

end
