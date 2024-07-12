function out = str(obj)
% Return string representation.

data = cell(1,3);
I = true(1,3);

% size
data(1) = compose('%dx%d',size(obj));

% nonzeros
data(2) = compose('%dnz',nnz(obj));

% degrees
d1 = obj.mindeg; d2 = obj.maxdeg;
if all(d2 == 0)
    % do not show degrees
    I(3) = false;
    % omit nonzeros for dense matrix
    I(2) = (nnz(obj) < numel(obj));
elseif d1 == d2
    % show single degree
    data(3) = compose('d=%d',d1);
else
    % show degree range
    data(3) = compose('d=%d:%d',d1,d2);
end

% join
out = join(data(I),',');

end
