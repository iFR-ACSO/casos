function [S,coeffs] = coeff_cat(obj,S2,coeff1,coeff2,dim)
% Concatenate polynomial coefficient matrices.

assert(any(dim==[1 2]),'Operating dimension must be either 1 or 2.')

% detect empty polynomials
if isempty(obj)
    % return second coefficient matrix
    S = S2;
    coeffs = coeff2;
    return

elseif isempty(S2)
    % return first coefficient matrix
    S = casos.Sparsity.create(obj);
    coeffs = coeff1;
    return
end

% else
sza = size(obj);
szb = size(S2);

assert(sza(3-dim) == szb(3-dim),'Dimensions of polynomials being concatenated are not consistent.')

% size of output
sz(dim) = sza(dim) + szb(dim); sz(3-dim) = sza(3-dim);

% combine variables
[indets,dga,dgb] = combineVar(obj.indets,S2.indets,obj.degmat,S2.degmat);

% extend coefficient matrices
% let a1,...,aL and b1,...,bM be the columns of A and B
[r1,c1] = get_triplet(obj.coeffs);
[r2,c2] = get_triplet(S2.coeffs);

switch (dim)
    case 1
        % vertical concatenation
        % vec([A;B]) = [(a1 b1) ... (aL bM)], assuming L=M

        % horizontal offset of i-th polynomial's 1st column
        offset1 =          szb(1)*floor(c1/sza(1));
        offset2 = sza(1) + sza(1)*floor(c2/szb(1));

    case 2
        % horizontal concatenation
        % vec([A B C]) = [a1...aL b1...bM c1...cN]

        % horizontal offset of i-th polynomial's coefficient matrix
        offset1 = 0;
        offset2 = sza(1).*sza(2);
end

% sparsity patterns of extended coefficient matrices
S1_coeffs = casadi.Sparsity.triplet(obj.nterm,prod(sz),r1,offset1+c1);
S2_coeffs = casadi.Sparsity.triplet(S2.nterm,prod(sz),r2,offset2+c2);

if isa(coeff1,'casadi.Sparsity')
    % concatenate sparsity
    coeffs = vertcat(S1_coeffs,S2_coeffs);

else
    % extend coefficient matrices
    cf1 = sparsity_cast(coeff1,S1_coeffs);
    cf2 = sparsity_cast(coeff2,S2_coeffs);
    % concatenate coefficients
    coeffs = vertcat(cf1,cf2);
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, [dga;dgb]);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,sz);

end
