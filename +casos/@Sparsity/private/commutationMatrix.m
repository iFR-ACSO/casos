function K = commutationMatrix(~, n)
% commutationMatrix: Construct commutation matrix K_{n,n}.
%   K = commutationMatrix(n) returns the sparse commutation matrix of size
%   n^2-by-n^2 such that:
%       K * vec(X) = vec(X.')

idx = 1:(n^2); 
[ii,jj] = ind2sub([n n],idx); 
idxT = sub2ind([n n],jj,ii); 
K = sparse(idxT,idx,1,n^2,n^2);

end