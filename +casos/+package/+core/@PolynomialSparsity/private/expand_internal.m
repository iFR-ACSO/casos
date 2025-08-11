function [cf1,cf2,degmat,indets] = expand_internal(obj,S2,cfa,cfb,omit_cf2)
% Internal function for expansion of polynomial coefficient matrices.

% combine variables
[indets,dga,dgb] = combineVar(obj.indets,S2.indets,obj.degmat,S2.degmat);

% make degree matrix unique
[S_cf1,degmat,I] = uniqueDeg(obj.coeffs, [dga;dgb]);

% expand coefficients
cf1 = sparsity_cast(cfa,S_cf1);

if nargin < 5 || ~omit_cf2
    % expand second coefficient matrix
    [id,ic] = I{:};
    [ii,jj] = get_triplet(S2.coeffs);
    S_cf2 = casadi.Sparsity.triplet(length(id),numel(S2),ic(obj.nterm+ii+1)-1,jj);
    cf2 = sparsity_cast(cfb,S_cf2);
else
    cf2 = [];
end

end
