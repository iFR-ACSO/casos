classdef QuasiconvBisection < casos.package.solvers.SosoptCommon
% Solve quasiconvex sum-of-squares problems via bisection.

properties (Access=private)
    sossolver;

    qc_sign;

    info = struct('iter',[]);
    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;

    log;
end

properties (Constant,Access=protected)
    bisect_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'conf_interval', 'Initial interval for the bisection.'
         'max_iter', 'Maximum number of bisections.'
         'sossol', 'The convex sum-of-squares solver to be used in the bisection.'
         'sossol_options', 'Options to be passed to the SOS solver.'
         'tolerance_abs', 'Absolute tolerance for stopping criterion.'
         'tolerance_rel', 'Relative tolerance for stopping criterion.'
         'verbose', 'Turn on/off iteration display.'}
    ];
end

properties (SetAccess=private)
    class_name = 'QuasiconvBisection';
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.QuasiconvBisection.bisect_options;
    end
end

methods (Access=protected)
    % internal evaluation
    argout = eval_on_basis(obj,argin);
end

methods
    function obj = QuasiconvBisection(name,sos,varargin)
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
            error('No objective given for quasiconvex problem.')
        else
            sos.f = casos.PS(sos.f);
    
            assert(isscalar(sos.f) && is_zerodegree(sos.f),'Objective must be scalar variable.')
            assert(is_symbolic(sos.f) || is_symbolic(-sos.f),'Objective must be symbolic variable.')
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
        if ~isfield(obj.opts,'conf_interval'), obj.opts.conf_interval = [-1e2 1e2]; end
        if ~isfield(obj.opts,'max_iter'), obj.opts.max_iter = 1000; end
        if ~isfield(obj.opts,'tolerance_abs'), obj.opts.tolerance_abs = 1e-3; end
        if ~isfield(obj.opts,'tolerance_rel'), obj.opts.tolerance_rel = 1e-3; end
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
