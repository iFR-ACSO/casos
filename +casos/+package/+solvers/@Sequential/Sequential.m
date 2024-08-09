classdef Sequential < casos.package.solvers.SosoptCommon
% Solve sum-of-squares problems via sequential SOS.

properties (Access=private)
    sossolver;

    % convex subproblem
    sizeHessian % needed for the initilization
    
    % second-order-correction
    solver_soc
    correctFun
    zeroDsoc
    
    % feasibility restoration phase
    solver_feas_res
    s0
    conVio_0
    feas_res_cost
    size_s
    size_x

    % BFGS
    BFGS_fun % function to efficiently evaluate BFGS

    % Filter
    Filter

    % parameterized projection (constraint violation)
    paraProjConVio
    xk1fun
   
    % cost function
    f
    nabla_xi_f
    nabla_xi_g 

    % constraint violation
    projConPara 
    pseudoProj


    % Langrangian
    nabla_xi_L_norm
    nabla_xi_L

    info = struct('iter',[]);
    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;

    log;
end

properties (Constant,Access=protected)
    sequential_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'max_iter', 'Maximum number of bisections.'
         'sossol', 'The convex sum-of-squares solver to be used in the subproblem.'
         'sossol_options', 'Options to be passed to the SOS solver.'
         'alpha_max', 'Maximum default step length'
         'alpha_min', 'Minimum default step length'
         'tau', 'Multiplier to adjust step length during linesearch'
         'soc_max_iter', 'Maximum number of sub-iterations in second-order correction'
         'optTol','Tolerance for optimality ||L_xi || <= optTol '
         'conVioTol','Tolerance for constraint violation || \theta(xi_k1) ||_2 <= conVioTol '
         'accTol','Tolerance for problem to be solved to acceptable level ( max(||L_xi ||, conVio) <= optTol )'
         'noAccIter','Number of iterations where the tolerance to acceptable level must be reached to leave optimization'
         's_phi', 'Cost function exponent for sufficient decrease condition (filter)'
         's_theta', 'Constraint violation exponent for sufficient decrease condition (filter)'
         'gamma_phi','Envelope parameter cost (filter)'
         'delta', 'Multiplier for sufficient decrease condition'
         'conViolCheck','Decide how to check for constraint violations.'
         'conVioSamp','Sampling points provided from user for pseudo-projection'
         'indeterminates','Vector of indeterminates'
         'verbose', ['Turn on/off iteration display.' ...
         'feasRes_actv_flag','flag to turn on/off feasibility restoration']}
    ];

    allow_eval_on_basis = true;
end

properties (SetAccess=private)
    class_name = 'Sequential';
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.Sequential.sequential_options;
    end
end

methods (Access=protected)
    % internal evaluation
    argout = eval_on_basis(obj,argin);
end

methods
    function obj = Sequential(name,sos,varargin)
        obj@casos.package.solvers.SosoptCommon(name,sos,varargin{:});
    
        % states
        if ~isfield(sos,'x')
            error('No decision variables given.')
        else
            sos.x = casos.PS(sos.x);
        end
        % parameter
        if ~isfield(sos,'p')
            % not parametrized
            sos.p = casos.PS;
        else
            sos.p = casos.PS(sos.p);
        end
        % objective
        if ~isfield(sos,'f')
            error('No objective given for sequential problem.')
        else
            sos.f = casos.PS(sos.f);
    
            assert(isscalar(sos.f) && is_zerodegree(sos.f),'Objective must be scalar variable.')
            % assert(is_symbolic(sos.f) || is_symbolic(-sos.f),'Objective must be symbolic variable.')
        end
        % constraints
        if ~isfield(sos,'g')
            error('No constraints given.')
        else
            sos.g = casos.PS(sos.g);
        end

        % default options
        if ~isfield(obj.opts,'sossol'),         obj.opts.sossol               = 'mosek'; end
        if ~isfield(obj.opts,'sossol_options'), obj.opts.sossol_options       = struct; end
        if ~isfield(obj.opts,'max_iter'),       obj.opts.max_iter             = 1000; end
        if ~isfield(obj.opts,'alpha_max'),      obj.opts.alpha_max            = 1; end
        if ~isfield(obj.opts,'alpha_min'),      obj.opts.alpha_min            = 1e-5; end
        if ~isfield(obj.opts,'tau'),            obj.opts.tau                  = 0.5; end
        if ~isfield(obj.opts,'soc_max_iter'),   obj.opts.soc_max_iter         = 5; end
        if ~isfield(obj.opts,'optTol'),         obj.opts.optTol               = 1e-4; end
        if ~isfield(obj.opts,'conVioTol'),      obj.opts.conVioTol            = 1e-2; end
        if ~isfield(obj.opts,'accTol'),         obj.opts.accTol               = 1e-2; end
        if ~isfield(obj.opts,'noAccIter'),      obj.opts.noAccIter            = 15; end
        if ~isfield(obj.opts,'s_phi'),          obj.opts.s_phi                = 2.3; end
        if ~isfield(obj.opts,'s_theta'),        obj.opts.s_theta              = 1.1; end
        if ~isfield(obj.opts,'gamma_phi'),      obj.opts.gamma_phi            = 1; end
        if ~isfield(obj.opts,'delta '),         obj.opts.delta                = 1; end

        if ~isfield(obj.opts,'conViolCheck'),     obj.opts.conViolCheck       = 'projection'; end
        if ~isfield(obj.opts,'conVioSamp'),     obj.opts.conVioSamp           = []; end
        if ~isfield(obj.opts,'indeterminates'), obj.opts.indeterminates       = []; end       
        if ~isfield(obj.opts,'feasRes_actv_flag'), obj.opts.feasRes_actv_flag = 1; end  


        % set up logger
        if ~isfield(obj.opts,'verbose') || ~obj.opts.verbose
            % no display
            obj.log = casos.package.Logger.Off;
        else
            % display debug messages
            obj.log = casos.package.Logger.Debug;
        end
    
        % build SOS problem
        buildproblem(obj,sos);
    end

    function s = get_stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end
end

end
