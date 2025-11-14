function K = blockCommutation(obj,Msdd)
% Mdd = [n1, n2, ..., nk] where each ni is a block size
% Build commutation matrices for each block

% Forces M to be a row vector
Msdd = Msdd(:).';   

% number of elements in each block 
Mq = Msdd.^2; 

% number of elements *before* block 
S = cumsum([0 Mq(1:end-1)]); 

% enumerate all elements 
idx = 1:sum(Mq); 

% compute index of block (1-based) 
J = sum(idx > S',1); 
% compute linear indices within each block 
idx0 = idx - S(J); 

% compute sub-indices for each block 
l = ceil(idx0 ./ Msdd(J)); % col 
k = idx0 - Msdd(J).*(l-1); % row 

% commute indices 
idx0T = (k-1).*Msdd(J) + l; 

% compute commutation of linear indices 
idxT = idx0T + S(J); 

% block commutation matrix  
K = sparse(idxT,idx,1);

end