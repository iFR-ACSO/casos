function V = sdp_vec(M, scale)
% Clarabel-style matrix vectorization for semidefinite cone embedding.
%
% This function takes a matrix
%
%   M = [ m1 ... ml ]
%
% where each column satisfies mi = Mi(:) for some k-by-k matrix Mi, and
% returns the matrix
%
%   V = [ v1 ... vl ]
%
% where each column is a k*(k+1)/2-by-1 vector that corresponds to stacking 
% the upper-triangular elements of Mi column-wise, where the off-diagonal 
% entries are scaled by sqrt(2).

k = sqrt(size(M,1));

assert(k == floor(k), 'Input must be a quadratic matrix.')

% user-defined scaling
if nargin < 2
    scale = sqrt(2);
end

% indices
[i, j] = ind2sub([k k], 1:k^2);
% transpose
idx = sub2ind([k k], j, i);
% ensure matrix is symmetric
M = (M(idx,:) + M)/2;

% select diagonal entries
Id = find(casadi.Sparsity.diag(k));
% upper-triangular entries
Iu = find(casadi.Sparsity.upper(k));

% scale upper-triangular off-diagonal entries
Msc = scale * M;
Msc(Id,:) = Msc(Id,:) + (1 - scale) * M(Id,:);

% stack elements column-wise
V = Msc(Iu,:); %#ok<FNDSB> % logical indexing not supported by CasADi

end