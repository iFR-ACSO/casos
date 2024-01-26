classdef Sequential < casos.package.solvers.SosoptCommon
% Solve nonlinear sum-of-squares problems via sequential sum-of-squares.

properties (Access=private)
    sossolver;
    lineSearch;

    info = struct('iter',[]);
    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
end


properties (Access=public)
   MeritFunEval;

end


properties (Constant,Access=protected)
    nlsos_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'max_iter', 'Maximum number of sequential sos.'
         'sossol', 'The convex sum-of-squares solver to be used in the bisection.'
         'sossol_options', 'Options to be passed to the SOS solver.'
         'tolerance_abs', 'Absolute tolerance for stopping criterion.'
         'tolerance_rel', 'Relative tolerance for stopping criterion.'}
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
        if ~isfield(obj.opts,'max_iter'), obj.opts.max_iter = 1000; end
        if ~isfield(obj.opts,'tolerance_abs'), obj.opts.tolerance_abs = 1e-3; end
        if ~isfield(obj.opts,'tolerance_rel'), obj.opts.tolerance_rel = 1e-3; end
    
       
        % store user parameter
        p0 = sos.p;
        
        % parameterize decision variables
        x0    = [];
        x1    = [];
        lam_x = [];
        lam_g = [];

        % for each decision variable generate a parameter (which is the current solution)
        % naming is not important since this is only for internal 
        for k = 1:length(sos.x)

            base = basis(sos.x(k));
           x0 = [x0; 
                         casos.PS.sym(['X0' num2str(k)],base)];

           x1 = [x1; 
                 casos.PS.sym(['X1' num2str(k)],base)];

           lam_x = [lam_x; 
                   casos.PS.sym(['lam_x' num2str(k)],base)];
        end


         for k = 1:length(sos.g)
           lam_g = [lam_g; 
                 casos.PS.sym(['lam_g' num2str(k)],basis(sos.g(k)))];
        end

        
        d = casos.PS.sym('d');
  
        ObjFun   = casos.Function('f',{sos.x}, {sos.f});
        ConFun   = casos.Function('f',{sos.x}, {sos.g});

        eta = 1e-15;


        %% setup line search solver
        % define SOS problem:
        %   min_d Psi(d) s.t. 0 \leq d \leq 1 
        sos_lineSearch = struct('x',d, ...
                                 'f',ObjFun((1-d)*x0 + d*x1) - dot(lam_x,x0) - dot(lam_g, ConFun((1-d)*x0 + d*x1)) - eta*d.^2, ...
                                 'g',[], ...
                                 'p',[x0; x1; lam_x; lam_g]);

        % constraint is scalar SOS cone
        opts = struct('Kc',struct('s',0),'Kx',struct('l',1));
        
        % solve by relaxation to SDP
        obj.lineSearch = casos.sossol('S','mosek',sos_lineSearch,opts);
       
         
        % Taylor Approximation constraints and evaluate at parameterized
        % solution
         sos.g = linearize(sos.g,sos.x,x0);

        % Taylor Approximation cost and evaluate at parameterized
        % solution
        sos.f = linearize(sos.f,sos.x,x0);
        
        % extend parameter vector
        sos.p = [p0; x0];

   
        
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
