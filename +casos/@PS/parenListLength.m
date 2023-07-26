function n = parenListLength(obj,indexOp,context)
% Overwriting matlab.mixin.indexing.RedefinesParen.parenListLength

import matlab.indexing.IndexingOperationType

if length(indexOp) < 2
    % indexing into polynomial matrix
    n = 1;
    return
end

% else
n = listLength(obj,indexOp(2:end),context);

end
