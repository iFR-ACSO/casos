classdef (Abstract) SequentialCommon < casos.package.solvers.SosoptCommon
    % Base class for sequential sum-of-squares algorithms.

    properties (Constant,Access=protected)
        sequential_options = [casos.package.solvers.SosoptCommon.sosopt_options
            {'sossol', 'The convex sum-of-squares solver to be used in the subproblem.'
            'sossol_options', 'Options to be passed to the SOS solver.'
            'tolerance_con', 'Absolute tolerance for stopping criterion of constraint violation.'
            'tolerance_opt', 'Absolute tolerance for stopping criterion of constraint violation.'
            'filter_struct', 'Structure containing parameter for filter linesearch.'
            'scale_BFGS0','Scaling parameter for initial BFGS matrix.'
            'Soc_is_enabled', 'Flag to turn on off the second-order-correction [default: true]'
            'Hessian_init','Method to initialize Hessian.'
            'hessian_approx','Hessian (Lagrangian) approximation method'
            'max_iter', 'Maximum number of iterations.'
            'almostOptCount','Number of iterations to check for almost optimal. '
            'feasibility_restoration','Control cost function of feasibility restoration.'
            'conVioCheck','Select method for constraint violation check'
            'userSample','If sampling is used, user can provide samples.'
            'verbose', 'Turn on/off iteration display.'}
            ];

        status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
    end

    properties (Access=protected)
        % low-level solvers
        solver_convex;
        solver_conVio;
        solver_soc;

        % Lagrangian and derivative
        eval_L;
        eval_dLdx;

        % linesearch
        eval_gradCost;
        eval_constraintSamples;

        % damped BFGS
        damped_BFGS;
        eval_s;
        eval_y;
        eval_r;

        % convergence check
        eval_cost;
        eval_constraint;
        eval_gradLag;
        eval_gradLag2;
        eval_Hessian;

        % parameters
        init_para;
        feasRes_para;
        display_para;

        % logging & info
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
        % prepare problem structures
        buildproblem(obj,nlsos);
        
        % iteration for overloading
        varargout = run_iteration(varargin);
        varargout = do_single_iteration(varargin);

        varargout = check_suffDecrease(varargin);

        % internal evaluation
        argout = eval_on_basis(obj,argin);
    end

    methods (Static)
        % iteration for overloading
        varargout = regularize_Hessian(varargin);
        varargout = hessian_regularization(varargin);
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

            % defaut parameter
            filter_struct = struct();

            filter_struct.alpha_max       = 1;
            filter_struct.s_theta         = 0.9;
            filter_struct.s_phi           = 2;
            filter_struct.gamma_theta     = 1e-5;
            filter_struct.gamma_phi       = 1e-5;
            filter_struct.gamma_alpha     = 0.05;
            filter_struct.delta           = 1;
            filter_struct.eta             = 1e-4;
            filter_struct.LangrangeFilter = 0;
            filter_struct.theta_min       = 1e-4;

            % make sure cholseky option is passed
            % sossol_options0.sdpsol_options.cholesky_method = 'num';

            % default options
            if ~isfield(obj.opts,'sossol'), obj.opts.sossol                                     = 'mosek'; end
            if ~isfield(obj.opts,'sossol_options'), obj.opts.sossol_options                     = struct; end
            if ~isfield(obj.opts,'tolerance_con'), obj.opts.tolerance_con                       = 1e-6; end
            if ~isfield(obj.opts,'tolerance_opt'), obj.opts.tolerance_opt                       = 1e-4; end
            if ~isfield(obj.opts,'scale_BFGS0'), obj.opts.scale_BFGS0                           = 1; end
            if ~isfield(obj.opts,'Hessian_init'),   obj.opts.Hessian_init                       = 'Analytical'; end
            if ~isfield(obj.opts,'hessian_approx'), obj.opts.hessian_approx                     = 'Regularization'; end
            if ~isfield(obj.opts,'feasibility_restoration'), obj.opts.feasibility_restoration   = 'Regularize'; end
            if ~isfield(obj.opts,'Soc_is_enabled'), obj.opts.Soc_is_enabled                     = true; end
            if ~isfield(obj.opts,'filter_struct'), obj.opts.filter_struct                       = filter_struct; end
            if ~isfield(obj.opts,'max_iter'), obj.opts.max_iter                                 = 100; end
            if ~isfield(obj.opts,'almostOptCount'), obj.opts.almostOptCount                     = 100; end
            if ~isfield(obj.opts,'conVioCheck'), obj.opts.conVioCheck                           = 'signed-distance'; end
            if ~isfield(obj.opts,'userSample'), obj.opts.userSample                             = []; end

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
