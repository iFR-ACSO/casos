function [cf1,cf2] = coeff_expand(obj,S2,cfa,cfb)
% Expand polynomial coefficient matrices to degrees.

% expand coefficients
[cf1,cf2] = expand_internal(obj,S2,cfa,cfb,(nargout < 2));

end
