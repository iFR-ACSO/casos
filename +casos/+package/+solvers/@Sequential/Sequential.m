classdef Sequential < casos.package.solvers.SosoptCommon
% Solve nonlinear sum-of-squares problems via sequential sum-of-squares.

properties (Access=private)
    sossolver;
    lineSearch;
        
    Merit;
    constraintFun
    cost_fun
    
    nabla_x_fun;
    nabla_lam_fun;

    info = struct('iter',[]);
    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
end


properties (Constant,Access=protected)
    nlsos_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'max_iter'      , 'Maximum number of sequential sos.'
         'sossol'        , 'The convex sum-of-squares solver to be used in the bisection.'
         'sossol_options', 'Options to be passed to the SOS solver.'
         'tolerance_abs' , 'Absolute tolerance for stopping criterion.'
         'tolerance_rel' , 'Relative tolerance for stopping criterion.'
         'verbose'       , 'Print current iteration to command window'
         'line_search'   , 'Select an algorithm to solve the linesearch problem.'}
    ];
end

properties (SetAccess=protected)
    class_name = 'Sequential';
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.Sequential.nlsos_options;
    end
end

methods
    out = call(obj,in);

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
            error('No objective given')
        else
            if ~isempty(sos.f)
                sos.f = casos.PS(sos.f);
            else
                sos.f = casos.PS(0);
            end
    
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
        if ~isfield(obj.opts,'tolerance_abs'), obj.opts.tolerance_abs = 1e-3; end
        if ~isfield(obj.opts,'tolerance_rel'), obj.opts.tolerance_rel = 1e-3; end
        if ~isfield(obj.opts,'verbose'), obj.opts.verbose = 0; end
         if ~isfield(obj.opts,'line_search'), obj.opts.line_search = 'bisection'; end
       
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
