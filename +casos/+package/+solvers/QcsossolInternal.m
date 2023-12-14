classdef QcsossolInternal < casos.package.functions.FunctionWrapper
% Internal interface for quasiconvex sum-of-squares problems.

methods
    function obj = QcsossolInternal(name,solver,varargin)
        % Create new sum-of-squares interface.
        switch (solver)
            case 'bisection'
                wrap = casos.package.solvers.QuasiconvBisection(name,varargin{:});
            otherwise
                error('No such quasiconvex solver "%s".',solver)
        end

        obj@casos.package.functions.FunctionWrapper(wrap);
    end
end

end
