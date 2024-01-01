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

elseif all(n==0)
    % power of zero is one
    b = casos.PS.ones(sz);
    return
end
% TODO: handle or escape for other simple cases, e.g., scalar, constant
% matrix, single term etc.?

% else
b = casos.PS;

nta = a.nterm; 
nva = a.nvars;
neb = prod(sz);
nen = numel(unique(n)); % count unique exponents


% reshape to output dimensions
deg = reshape(repmat(n,sz./szn),1,prod(sz));
cfa = reshape(repmat(a.coeffs,sz./sza),nta,prod(sz));
dga = a.degmat;

% check number of coefficients per matrix component
idx = find(sparsity(cfa));
[~,ja] = ind2sub(size(cfa),idx);

if nva == 0 ... % base polynomial is a constant
    || (nta == 1 && nen == 1) % a^n = c*(x^d)^n
    b.coeffs = cfa.^deg;
    b.degmat = n.*dga;
    b.indets = a.indets;
    b.matdim = sz;
    return

elseif issorted(ja,'strictascend')
    % base polynomial is matrix of monomials
    % match exponents to unique degrees
    [dd,~,Ideg] = unique(deg);
    % repeat coefficients and degrees to match exponents
    cfa = repmat(cfa,nen,1);
    dga = repmat(dga,nen,1);
    % match coefficients to corresponding monomials
    [ia,ja] = get_triplet(sparsity(cfa)); % CasADi interface has 0-index
    is_deg = ceil((ia(:)+1)./nva) == Ideg(ja(:)+1);
    % select matching coefficients
    S = casadi.Sparsity.triplet(size(cfa,1),size(cfa,2),ia(is_deg),ja(is_deg));
    coeffs = project(cfa,S).^repmat(deg,nva*nen,1);
    % multiply degree matrix with (unique) exponents
    degmat = dga.*reshape(repmat(dd,nva,1),nen*nva,1);

else
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
C(1) = {casadi.SX.ones(1,neb)}; D(1) = {sparse(1,nva)};
C(2) = {cfa};                   D(2) = {dga};

[C,D] = powers(C,D,max(n),nta);

if nen > 1
    % either matrix.^matrix or vector.^vector
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
        Cd{i} = project(cfi, S);
    end

    coeffs = vertcat(Cd{:});
    degmat = vertcat(D{1+deg});
else
    % matrix.^scalar
    coeffs = C{end};
    degmat = D{end};
end

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