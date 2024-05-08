function [S,cf1,cf2] = coeff_expand(S1,S2,cfa,cfb)
% Expand polynomial coefficient matrices to degrees.

S = casos.Sparsity;

% combine variables
[indets,dga,dgb] = combineVar(S1.indets,S2.indets,S1.degmat,S2.degmat);

% all-sparse zero matrix
zfb = getZero(cfa,cfb);

% make degree matrix unique
[cf1,degmat,I] = uniqueDeg([cfa;zfb], [dga;dgb]);

if nargout > 2
    % expand second coefficient matrix
    zfa = getZero(cfb,cfa);
    cf2 = uniqueDeg([zfa;cfb], I);
end

% new sparsity pattern
S.degmat = degmat;
S.indets = indets;
S.matdim = size(S1);
% store coefficients
S = set_coefficients(S,cf1);

end

function z = getZero(c1,c2)
% Return an all-sparse zero matrix.

    if isa(c1,'casadi.Sparsity')
        % all-sparse pattern
        z = casadi.Sparsity(size(c2));
    
    else
        % all-sparse matrix
        z = sparse(size(c2,1),size(c2,2));
    end
end