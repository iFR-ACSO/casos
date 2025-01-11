classdef (Abstract) SequentialCommon < casos.package.solvers.SosoptCommon
% Base class for sequential sum-of-squares algorithms.

properties (Constant,Access=protected)
    sequential_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'sossol', 'The convex sum-of-squares solver to be used in the bisection.'
         'sossol_options', 'Options to be passed to the SOS solver.'
         'tolerance_abs', 'Absolute tolerance for stopping criterion.'
         'tolerance_rel', 'Relative tolerance for stopping criterion.'
         'max_iter', 'Maximum number of iterations.'
         'verbose', 'Turn on/off iteration display.'}
    ];

    allow_eval_on_basis = true;


    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
end

% public properties to allow e.g. feasibility restoration accesses to it
properties

        solver_conVio

        % functions to be evaluated (convergence check)
        eval_cost
        hess_fun
        % linesearch
        eval_gradCost
        
        % Langrangian and derivative
        L
        dLdx 

        init_para

        FeasRes_para
        
        sparsity_pat_para

end

properties (Access=protected)
    % low-level solver
    solver_convex;
    solver_soc

    % damped BFGS
    damped_BFGS
    SR1
    eval_s
    eval_y
    eval_r

    % functions to be evaluated (convergence check)
    % eval_cost
    eval_gradLang
    eval_gradLang2


    % parameter/data for display output
    display_para
end

properties (Access=protected)
    log;
    info = struct('iter',[]);
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SequentialCommon.sequential_options;
    end
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
    function obj = SequentialCommon(name,nlsos,varargin)
        obj@casos.package.solvers.SosoptCommon(name,nlsos,varargin{:});

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
