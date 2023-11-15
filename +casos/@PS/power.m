function b = power(a,n)
% Element-wise powers.

a = casos.PS(a);

% input dimensions
sza = size(a);
szn = size(n);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% compare dimensions
I = (sza == szn);
% find zero dimension
I0 = (sza == 0) | (szn == 0);
% find one dimension
I1 = (sza == 1) | (szn == 1);

% dimensions are compatible if equal or one factor is row/column
assert(all(I | I1), errsz, size2str(sza), size2str(szn))

% dimensions of element-wise product
sz = max(sza,szn);

% handle simple case(s) for speed up
if isempty(a) || isempty(n)
    % element-wise power with empty polynomial/exponent is empty
    sz(I0) = 0;

    b = casos.PS.zeros(sz);
    return
end
% TODO: handle or escape for other simple cases, e.g., scalar, constant
% matrix, single term etc.?

% else
b = casos.PS;

nta = a.nterm; 
nva = a.nvars;
nea = numel(a);
nen = numel(n);


% reshape to output dimensions
cfa = reshape(repmat(a.coeffs,sz./sza),nta,prod(sz));
dga = a.degmat;

% if nen > 1
%     % modify coefficients so every matrix component has a separate row
%     idx = find(sparsity(cfa));
%     [ia,ja] = ind2sub(size(cfa),idx);
%     % repeat column vectors for each component
%     S = sparsity(casadi.SX(sparse(ia + (ja-1)*nta,ja,1,nta*nea,nea)));
%     cfa = casadi.SX(S, cfa(idx));
% 
%     % repeat degree matrix to match
%     dga = repmat(a.degmat,nea,1);
% end

% compute powers
C = cell(1,max(n)+1);           D = cell(1,max(n)+1);
C(1) = {casadi.SX.ones(1,nea)}; D(1) = {sparse(1,nva)};
C(2) = {cfa};                   D(2) = {dga};

[C,D] = powers(C,D,max(n),nta);

if nen > 1
    % either matrix.^matrix or vector.^vector
    deg = repmat(n,sz./szn);
    % prepare components
    Cd = cell(size(deg));
    % iterate over powers
    for i=1:numel(deg)
        cfi = C{deg(i)+1};
        [mi,ni] = size(cfi);
        % isolate i-th column of coefficients
        idx = find(sparsity(cfi));
        [ii,ji] = ind2sub(size(cfi),idx);
        I = (ji == i);
        S = sparsity(casadi.SX(sparse(ii(I),ji(I),1,mi,ni)));
        Cd{i} = casadi.SX(S, cfi(:,i));
    end

    coeffs = vertcat(Cd{:});
    degmat = vertcat(D{1+deg});
else
    % matrix.^scalar
    coeffs = C{n+1};
    degmat = D{n+1};
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,a.indets);

% new polynomial
b.coeffs = coeffs;
b.degmat = degmat;
b.indets = indets;
b.matdim = sz;

end

function [C,D] = powers(C,D,l,nt)
% Compute coefficient and degree matrices up to power l.

if ~isempty(C{l+1})
    % nothing to do
    return

elseif rem(l,2) > 0
    % odd number
    [C,D] = powers(C,D,l-1,nt);

    nta = nt^(l-1);
    ntb = nt;
    
    coeffs = kron(C{l},ones(ntb,1)) .* kron(ones(nta,1),C{2});
    degmat = kron(D{l},ones(ntb,1)) +  kron(ones(nta,1),D{2});

else
    [C,D] = powers(C,D,l/2,nt);
    [C,D] = powers(C,D,l-1,nt);

    nta = nt^(l/2);
    ntb = nt^(l/2);
    
    coeffs = kron(C{1+l/2},ones(ntb,1)) .* kron(ones(nta,1),C{1+l/2});
    degmat = kron(D{1+l/2},ones(ntb,1)) +  kron(ones(nta,1),D{1+l/2});
end

C(l+1) = {coeffs};
D(l+1) = {degmat};

end