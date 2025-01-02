function [M,pts_i,pts_o] = op2interpol(op,pts_i,pts_o)
% Return a matrix of coordinates for given input and output interpolations.
% If no interpolations are given, the interpolation basis are computed 
% based on the operator's sparsity patterns.

if nargin < 2
    % compute interpolation basis
    [pts_i,Vi] = interpolation(op.sparsity_in);
    [pts_o,Vo] = interpolation(op.sparsity_out);

else
    % compute Vandermonde matrices
    Vi = vandermonde(op.sparsity_in, pts_i);
    Vo = vandermonde(op.sparsity_out, pts_o);
end

% compute linear operator from input interpolation to output interpolation
M = Vo*op.matrix/Vi;

end
