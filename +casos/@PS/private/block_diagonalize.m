function [K,Lz, nz] = block_diagonalize(Pdegmat, degmat, N)
% search for simple symmetries in the polynomial, and block-diagonalize 
% the problem

% obtain degree matrices of polynomials and monomial basis
P = Pdegmat';

% check if degmat is empty
if isempty(degmat)
    return 
end

% obtain size of indeterminates and matrix R
allsym = 2^N;
R = dec2bin(0:allsym-1) - '0';
idx = 1:allsym;

% find symmetries
signR = arrayfun(@(i) ~any(mod(R(i,:)*P, 2)), idx, 'UniformOutput', false);
signR = vertcat(signR{:});

% remove non-sign symmetric rows
R(~signR,:)=[];
W = R*degmat';
W = mod(W,2);
[~,~,IC]=unique(W','rows');
[gc,idx] = groupcounts(IC);
K = gc';

% Obtain a Lz vector that chooses the elements from the basis
Lz = arrayfun(@(i) IC==i, idx, 'UniformOutput', false);
Lz = horzcat(Lz{:})';
nz = length(idx);
end