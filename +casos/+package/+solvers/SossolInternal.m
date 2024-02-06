classdef SossolInternal < casos.package.functions.FunctionWrapper
% Internal interface for convex sum-of-squares problems.

methods
    function obj = SossolInternal(varargin)
        % Create new sum-of-squares interface.
        wrap = casos.package.solvers.SossdpRelaxation(varargin{:});

        obj@casos.package.functions.FunctionWrapper(wrap);
    end
end

end
