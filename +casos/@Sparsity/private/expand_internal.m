function [cf1,cf2,degmat,indets] = expand_internal(S1,S2,cfa,cfb,omit_cf2)
% Internal function for expansion of polynomial coefficient matrices.

% combine variables
[indets,dga,dgb] = combineVar(S1.indets,S2.indets,S1.degmat,S2.degmat);

% all-sparse zero matrix
zfb = getZero(cfa,cfb);

% make degree matrix unique
[cf1,degmat,I] = uniqueDeg([cfa;zfb], [dga;dgb]);

if nargin < 5 || ~omit_cf2
    % expand second coefficient matrix
    zfa = getZero(cfb,cfa);
    cf2 = uniqueDeg([zfa;cfb], I);
else
    cf2 = [];
end

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
