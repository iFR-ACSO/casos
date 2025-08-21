function [S,I] = op_subsref(obj,ii)
% Subreference into nonzero coefficients.

[S,~,I1,I2] = coeff_subsref(obj,obj.sparsity_M,ii,[]);

% get indices of referenced nonzeros
[ic,jc] = get_triplet(obj.sparsity_M);
% find referenced coefficients
I = find(ismember(I1,ic+1) & ismember(I2,jc+1));    % CasADi has 0-index

end
