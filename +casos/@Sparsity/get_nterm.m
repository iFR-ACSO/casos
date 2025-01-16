function nt = get_nterm(obj,I)
% Get number of terms.

[~,L] = get_degmat(obj);

if nargin < 2
    % return number of terms per element
    nt = sum(L,2);
else
    % return number of terms in selection
    nt = nnz(any(L(I,:),1));
end

end
