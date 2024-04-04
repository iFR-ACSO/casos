classdef (Abstract) PolynomialInterface
% Base class for polynomial-like objects.

methods (Abstract)
    % check if object represents indeterminate variables
    tf = is_indet(obj);
end

end
