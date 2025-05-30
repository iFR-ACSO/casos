function M = sdp_mat(obj,V,scale)
% SCS-style matrix de-vectorization for semidefinite cone embedding.
%
% This function takes a matrix
% 
%   V = [ v1 ... vl ]
% 
% where each column is a k*(k+1)/2-by-1 vector that corresponds to stacking 
% the (lower or upper) triangular elements of a k-by-k symmetric matrix Mi 
% column-wise, where the off-diagonal entries are scaled by sqrt(2), and 
% returns the matrix
%
%   M = [ m1 ... ml ]
%
% where each column satisfies mi = Mi(:).

k = (sqrt(8*size(V,1) + 1) - 1)/2;

assert(k == floor(k), 'Input must be a compatible vector.')

% user-defined scaling
if nargin < 3
    scale = sqrt(2);
end

% select triangular elements
s = obj.sparsity_triangular(k);
S = repmat(reshape(s,k^2,1),1,size(V,2));

% assign elements column-wise
Msc = feval(class(V),S,V(:));

% select diagonal entries
Id = find(casadi.Sparsity.diag(k));

% indices
[i,j] = ind2sub([k k],1:k^2);
% transpose
idx = sub2ind([k k],j,i);

% mirror off-diagonal triangular entries,
% note: this ensures output is symmetric
% de-scale off-diagonal entries
M = scale\(Msc(idx,:) + Msc);
M(Id,:) = M(Id,:) + (1 - 2/scale)*Msc(Id,:);

end
