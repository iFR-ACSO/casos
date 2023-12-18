function V = sdp_vec(M)
% MOSEK/SCS-style matrix vectorization for semidefinite cone embedding.
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
% the lower-triangular elements of Mi column-wise, where the off-diagonal 
% entries are scaled by sqrt(2).

k = sqrt(size(M,1));

assert(k == floor(k), 'Input must be a quadratic matrix.')

% indices
[i,j] = ind2sub([k k],1:k^2);
% transpose
idx = sub2ind([k k],j,i);
% ensure matrix is symmetric
M = (M(idx,:) + M)/2;

% select diagonal entries
Id = find(casadi.Sparsity.diag(k));
% lower-triangular entries
Il = find(casadi.Sparsity.lower(k));

% scale lower-triangular off-diagonal entries, 
Msc = sqrt(2)*M;
Msc(Id,:) = Msc(Id,:) + (1-sqrt(2))*M(Id,:);

% stack elements column-wise
V = Msc(Il,:);

end
