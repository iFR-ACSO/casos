classdef SimpleSequential < casos.package.solvers.SequentialCommon
% A simple sequential sum-of-squares algorithm.

properties (SetAccess=private)
    class_name = 'SimpleSequential';
end

properties(Access=protected)
    feas_res_solver
end

methods (Access=protected)
    % iteration for overloading
    varargout = run_iteration(varargin);
    varargout = do_single_iteration(varargin);

    % internal evaluation
    argout = eval_on_basis(obj,argin);
end

methods
    % Constructor
    function obj = SimpleSequential(name,nlsos,varargin)
        obj@casos.package.solvers.SequentialCommon(name,nlsos,varargin{:});

        % states
        if ~isfield(nlsos,'x')
            error('No decision variables given.')
        else
            nlsos.x = casos.PS(nlsos.x);
        end
        % parameter
        if ~isfield(nlsos,'p')
            % not parametrized
            nlsos.p = casos.PS;
        else
            nlsos.p = casos.PS(nlsos.p);
        end
        % objective
        if ~isfield(nlsos,'f')
            nlsos.f = casos.PS(0);
        else
            nlsos.f = casos.PS(nlsos.f);
        end
        % constraints
        if ~isfield(nlsos,'g')
            nlsos.g = casos.PS;
        else
            nlsos.g = casos.PS(nlsos.g);
        end

        % default options
        if ~isfield(obj.opts,'sossol'), obj.opts.sossol = 'mosek'; end
        if ~isfield(obj.opts,'sossol_options'), obj.opts.sossol_options = struct; end
        if ~isfield(obj.opts,'tolerance_abs'), obj.opts.tolerance_abs = 1e-3; end
        if ~isfield(obj.opts,'tolerance_rel'), obj.opts.tolerance_rel = 1e-3; end
        if ~isfield(obj.opts,'max_iter'), obj.opts.max_iter = 10; end
        % set up logger
        if ~isfield(obj.opts,'verbose') || ~obj.opts.verbose
            % no display
            obj.log = casos.package.Logger.Off;
        else
            % display debug messages
            obj.log = casos.package.Logger.Debug;
        end

        % Initialize solvers
        buildproblem(obj,nlsos);
    end

    function s = get_stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end

end
end




