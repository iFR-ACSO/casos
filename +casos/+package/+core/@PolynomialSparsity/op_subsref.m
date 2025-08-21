function [S,I] = op_subsref(obj,ii)
% Subreference into nonzero coefficients.

S = coeff_subsref(obj,obj.coeffs,ii,[length(ii) 1]);

% get indices of referenced nonzeros
[~,jc] = get_triplet(obj.coeffs);
% find referenced coefficients
I = find(ismember(ii,jc+1));    % CasADi has 0-index

end
