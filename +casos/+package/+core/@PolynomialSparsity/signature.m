function dim = signature(obj,never_empty)
% Return a signature representation of sparsity pattern.

[n,m] = size(obj);
d1 = obj.mindeg;
d2 = obj.maxdeg;

if n == 1 && m == 1 && (nargin < 2 || ~never_empty || d2 > 0)
    % don't show dimensions
    dim = {};

elseif m == 1
    % show number of rows
    dim = compose('%d',n);

elseif n > 0 || m > 0
    % show size
    dim = compose('%dx%d',n,m);

else
    % show empty dimensions
    dim = {''};
end

if all(d2 == 0)
    % nothing to do
    return

elseif ~isempty(dim)
    % add degrees
    dim{end+1} = ',';
end

if d1 == d2
    % show single degree
    dim{end+1} = sprintf('d=%d',d1);
else
    % show degree range
    dim{end+1} = sprintf('d=%d:%d',d1,d2);
end

end
