% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs] = coeff_power(S,coeffs,deg)
% Compute element-wise power of coefficient matrix.

n = unique(deg); % count unique exponents

nta = S.nterm; 
nva = S.nvars;
neb = S.numel;
nen = numel(n);

% reshape to dimensions of coefficient matrix
deg = reshape(deg,1,neb);

if nva == 0 ... % base polynomial is a constant
    || (nta == 1 && nen == 1) % a^n = c*(x^d)^n
    % find zero degrees
    in = find(deg == 0); % 1-index

    if ~isa(coeffs,'casadi.Sparsity')
        coeffs = coeffs.^deg;
        % set zero powers to one
        if ~isempty(in), coeffs(:,in) = 1; end
    else
        % copy pattern
        [ii,jj] = get_triplet(coeffs);
        % add zero powers (always dense)
        coeffs = casadi.Sparsity.triplet(nta,neb,[ii zeros(size(in))],[jj in-1]); % 0-index
    end
    if nva > 0
        S.degmat = n.*S.degmat;
    end
    % store coefficients
    S = set_coefficients(S,coeffs);
    return

elseif is_monom(S)
    % base polynomial is matrix of monomials
    % find zero degrees
    in = find(deg == 0); % 1-index
    
    % reorder coefficients to match powers
    [ii,jj] = get_triplet(S.coeffs); % 0-index
    % add zero powers (always dense)
    S_coeffs = casadi.Sparsity.triplet(neb,neb,[jj in-1],[jj in-1]); % 0-index

    if ~isa(coeffs,'casadi.Sparsity')
        % take power of nonzero coefficients
        cfa = sum1(coeffs).^deg;
        % set zero powers to one
        if ~isempty(in), cfa(in) = 1; end
        % cast onto reordered sparsity pattern
        coeffs = sparsity_cast(cfa,S_coeffs);
    else
        % use coefficient sparsity pattern
        coeffs = S_coeffs;
    end

    % reorder degrees to match powers
    dga = sparse(neb,nva);
    dga(jj+1,:) = S.degmat(ii+1,:); % 1-index
    % multiply degree matrix with exponents
    degmat = dga.*reshape(deg,neb,1);

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
    D(2) = {S.degmat};

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
[coeffs,S.degmat,S.indets] = removeZero(coeffs,degmat,S.indets);

% store coefficients
S = set_coefficients(S,coeffs);

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
