classdef Sequential < casos.package.solvers.SosoptCommon
% Solve nonlinear sum-of-squares problems via sequential sum-of-squares.

properties (Access=private)
    sossolver;
    lineSearch;
    projCon;
    idxNonlinCon;
    solverSOC
    socConstSoc
        
    Merit;
    constraintFun;
    cost_fun;
    conFunRed;
    cost
    delta_dual
    delta_search
    FeasResPhase 
    projConPara;
    newIterate
    
    nabla_f 
    nabla_x_fun;
    nabla_lam_fun;
    langrangeLinear
    
    samplingPoint_proj

    % helper functions
    plusFun;
    vertcatFun;
    updateLineSearch
    deltaOptVar;
    norm2FunOptVar
    norm2FunVio
    BFGS_fun
    base_x
    

    dLdx
    log;
    sizeHess;

    info = struct('iter',[]);
    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
end


properties (Constant,Access=protected)
    nlsos_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'max_iter'      , 'Maximum number of sequential sos.'
         'sossol'        , 'The convex sum-of-squares solver to be used in the sequential SOS.'
         'sossol_options', 'Options to be passed to the SOS solver.'
         'tolerance_abs' , 'Absolute tolerance for stopping criterion.'
         'tolerance_rel' , 'Relative tolerance for stopping criterion.'
         'verbose'       , 'Turn on/off iteration display.'
         'optimality_tol', 'Optimality tolerance'
         'dual_tol', 'Dual variable tolerance'
         'conVio_tol', 'Constraint violation tolerance'
         'globalization'   , 'Select an algorithm for globalization strategy.'
         'Sequential_Algorithm','Selecte between SLP and SQP like algorithms.'
         'Sampling_Point', 'Sampling points to check constraint violation'
         'FeasRes', 'Dummy Parameter to initialize and use FeasRes'}
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
        if ~isfield(obj.opts,'max_iter'), obj.opts.max_iter = 100; end
        if ~isfield(obj.opts,'optimality_tol'), obj.opts.optimality_tol = 1e-3; end
        if ~isfield(obj.opts,'conVio_tol'), obj.opts.conVio_tol = 1e-5; end
        if ~isfield(obj.opts,'dual_tol'), obj.opts.dual_tol = 1e-5; end
        if ~isfield(obj.opts,'Sequential_Algorithm'), obj.opts.Sequential_Algorithm = 'SLP'; end
        if ~isfield(obj.opts,'Sampling_Point'), obj.opts.Sampling_Point = zeros(2,1000); end
        if ~isfield(obj.opts,'FeasRes'), obj.opts.FeasRes = 0; end
        if ~isfield(obj.opts,'line_search'), obj.opts.globalization = 'filter_linesearch'; end % filter_linesearch_simple
       
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