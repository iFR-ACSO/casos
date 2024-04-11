classdef (Abstract) PolynomialInterface < casos.package.core.Printable
% Base class for polynomial-like objects.

methods (Access=protected)
    %% Utilities
    function str = size2str(obj)
        % Convert array size to string.
        str = sprintf('%dx%d', size(obj,1), size(obj,2));
    end
end

end
