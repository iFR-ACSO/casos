function [S,coeffs] = coeff_power(obj,coeffs,deg)
% Compute element-wise power of coefficient matrix.

n = unique(deg); % count unique exponents

nta = obj.nterm; 
nva = obj.nvars;
neb = obj.numel;
nen = numel(n);

% reshape to dimensions of coefficient matrix
deg = reshape(deg,1,neb);

if nva == 0 ... % base polynomial is a constant
    || (nta == 1 && nen == 1) % a^n = c*(x^d)^n
    if ~isa(coeffs,'casadi.Sparsity')
        coeffs = coeffs.^deg;
    else
        % copy pattern
        coeffs = casadi.Sparsity(coeffs);
    end
    % multiply degrees
    degmat = n.*obj.degmat;

    % new sparsity pattern
    S = new_from_coefficients(coeffs,degmat,obj.indets,size(obj));
    return

elseif is_monom(obj)
    % base polynomial is matrix of monomials
    % match exponents to unique degrees
    [dd,~,Ideg] = unique(deg);
    % repeat coefficient sparsity to match exponents
    S_cfa  = repmat(obj.coeffs,nen,1);
    % match coefficients to corresponding monomials
    [ia,ja] = get_triplet(S_cfa); % CasADi interface has 0-index
    is_deg = ceil((ia(:)+1)./nta) == Ideg(ja(:)+1);
    % select matching coefficients
    S_cfb = casadi.Sparsity.triplet(size(S_cfa,1),size(S_cfa,2),ia(is_deg),ja(is_deg));
    if ~isa(coeffs,'casadi.Sparsity')
        % repeat coefficients to match exponents
        cfa = repmat(coeffs,nen,1);
        % project onto new sparsity pattern
        coeffs = project(cfa,S_cfb).^repmat(deg,nta*nen,1);
    else
        % use new sparsity pattern
        coeffs = S_cfb;
    end
    % repeat degrees to match exponents
    dga = repmat(obj.degmat,nen,1);
    % multiply degree matrix with (unique) exponents
    degmat = dga.*reshape(repmat(dd,nta,1),nen*nta,1);

else
    % compute all (element-wise) powers
    n_max = max(n);

    % prepare powers of coefficients
    C = cell(1,n_max+1);
    if ~isa(coeffs,'casadi.Sparsity')
        C(1) = {coeffs.ones(1,neb)};
    else
        C(1) = {casadi.Sparsity.dense(1,neb)};
    end
    C(2) = {coeffs};

    % prepare powers of degree matrix
    D = cell(1,n_max);
    D(1) = {sparse(1,nva)};
    D(2) = {obj.degmat};

    % compute powers up to maximal degree
    [C,D] = powers(C,D,n_max,nta);

    if nen > 1
        % either matrix.^matrix or vector.^vector
        % prepare components
        Cd = cell(size(deg));
        % iterate over powers
        for i=1:numel(deg)
            cfi = C{deg(i)+1};
            [mi,ni] = size(cfi);
            % isolate i-th column of coefficients
            if ~isa(cfi,'casadi.Sparsity')
                % idx = find(sparsity(cfi));
                % [ii,ji] = ind2sub(size(cfi),idx);
                [ii,ji] = get_triplet(sparsity(cfi));
            else
                [ii,ji] = get_triplet(cfi);
            end
            I = (ji == i-1);
            Si = casadi.Sparsity.triplet(mi,ni,ii(I),ji(I));
            if ~isa(cfi,'casadi.Sparsity')
                Cd{i} = project(cfi, Si);
            else
                Cd{i} = intersect(cfi, Si);
            end
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
[coeffs,degmat,indets] = removeZero(coeffs,degmat,obj.indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,size(obj));

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
        
        coeffs = times_coeffs(kron(C{l},ones(ntb,1)), kron(ones(nta,1),C{2}));
        degmat = times_degmat(kron(D{l},ones(ntb,1)), kron(ones(nta,1),D{2}));
    
    else
        [C,D] = powers(C,D,l/2,nt);
        [C,D] = powers(C,D,l-1,nt);
    
        nta = nt^(l/2);
        ntb = nt^(l/2);
        
        coeffs = times_coeffs(kron(C{1+l/2},ones(ntb,1)), kron(ones(nta,1),C{1+l/2}));
        degmat = times_degmat(kron(D{1+l/2},ones(ntb,1)), kron(ones(nta,1),D{1+l/2}));
    end
    
    C(l+1) = {coeffs};
    D(l+1) = {degmat};
end
