classdef Sequential < casos.package.solvers.SosoptCommon
% Solve sum-of-squares problems via sequential SOS.

properties (Access=private)
    sossolver;

    % convex subproblem
    sizeHessian % needed for the initilization
    
    % second-order-correction
    solver_soc

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

    % constraint violation
    projConPara 

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
         'verbose', 'Turn on/off iteration display.'}
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
        if ~isfield(obj.opts,'sossol'), obj.opts.sossol = 'sedumi'; end
        if ~isfield(obj.opts,'sossol_options'), obj.opts.sossol_options = struct; end
        if ~isfield(obj.opts,'max_iter'), obj.opts.max_iter = 1000; end


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
