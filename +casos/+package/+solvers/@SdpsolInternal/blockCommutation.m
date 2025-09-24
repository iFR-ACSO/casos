function K = blockCommutation(~,Msdd)
% Mdd = [n1, n2, ..., nk] where each ni is a block size
% Build commutation matrices for each block
Kblocks = cellfun(@commutationMatrix, num2cell(Msdd), 'UniformOutput', false);

% Assemble block diagonal matrix
K = blkdiag(Kblocks{:});

% Ensure sparse format
K = sparse(K);
end