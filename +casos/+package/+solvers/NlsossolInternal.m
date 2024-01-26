classdef NlsossolInternal < casos.package.functions.FunctionWrapper
% Internal interface for quasiconvex sum-of-squares problems.

methods
    function obj = NlsossolInternal(name,solver,varargin)
        % Create new sum-of-squares interface.
        switch (solver)
            case 'sequential'
                wrap = casos.package.solvers.Sequential(name,varargin{:});
            otherwise
                error('No such nonlinear solver "%s".',solver)
        end

        obj@casos.package.functions.FunctionWrapper(wrap);
    end
end

end
