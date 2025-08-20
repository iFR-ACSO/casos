classdef FeasibilityRestoration < casos.package.solvers.SequentialCommon
    % A simple sequential sum-of-squares algorithm.

    properties (SetAccess=private)
        class_name = 'FeasibilityRestoration';
    end

    properties(Access=protected)
        filter_feas
    end

    methods (Access=protected)
        % iteration for overloading
        varargout = run_iteration(varargin);
        varargout = do_single_iteration(varargin);

        % internal evaluation
        argout = eval_on_basis(obj,argin);

    end

    methods
        varargout = eval_extended(varargin);

    end



end
