function S = sub(obj,ii,S0)
% Subreference polynomial sparsity pattern.

if ~isa(S0,'casos.Sparsity') && ~isa(S0,'casadi.Sparsity')
    % subindices given
    [R,C] = ndgrid(ii,S0);
    ii = sub2ind(size(obj),R(:),C(:));
    sz = size(R);

else
    assert(length(ii) == numel(S0), 'Not supported.')
    sz = size(S0);
end

% subreference coefficient pattern
S = coeff_subsref(obj,obj.coeffs,ii,sz);

end
