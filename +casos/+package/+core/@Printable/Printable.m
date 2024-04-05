classdef (Abstract) Printable
% Base class for display.

methods (Abstract)
    % Return string representation(s) as cell array.
    out = str(obj);
end

methods
    function disp(obj)
        % Print string representation to command line.
        disp(cell2mat(str(obj)));
    end
end

methods (Access=protected)
    disp_matrix(obj,varargin);
end

end
