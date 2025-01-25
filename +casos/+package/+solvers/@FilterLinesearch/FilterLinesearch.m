classdef FilterLinesearch < casos.package.solvers.SequentialCommon
% A simple sequential sum-of-squares algorithm.

properties (SetAccess=private)
    class_name = 'FilterLinesearch';
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
    function obj = FilterLinesearch(name,nlsos,varargin)
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

        %% set up feasibility restoration phase --> do we need all or only nonlinear
        
        base_g = sparsity(nlsos.g);
        base_x = sparsity(nlsos.x);
        
        r  = casos.PS.sym('r',length(nlsos.g));
        s0 = casos.PD(base_g,ones(base_g.nnz,1));
        
        nlsos_feas.x = [r;nlsos.x];
        
        nlsos_feas.g = nlsos.g + r.*s0;
        
        x_R   = casos.PS.sym('x_R',base_x);
        nlsos_feas.f = sum(r) ; %+ 0.1/2*dot(nlsos.x-x_R,nlsos.x-x_R);
        
        nlsos_feas.p = [nlsos.p; x_R];
        
        obj.FeasRes_para.n_r = length(nlsos.g);
        
        
        sosopt               = obj.opts.sossol_options;
        sosopt.sossol        = obj.opts.sossol;
        sosopt.Kx.lin        = length(nlsos_feas.x);
        sosopt.Kc.sos        = length(nlsos_feas.g);
        sosopt.error_on_fail = false;
        sosopt.verbose       = 1;
        sosopt.max_iter      = 100;
        obj.feas_res_solver  =  casos.package.solvers.FeasibilityRestoration('feasRes',nlsos_feas,sosopt);
        
       % total build time for both actual problem and feasibility
       % restoration
       obj.display_para.solver_build_time      = obj.display_para.solver_build_time +obj.feas_res_solver.display_para.solver_build_time;
    end

    function s = get_stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end

end
end




