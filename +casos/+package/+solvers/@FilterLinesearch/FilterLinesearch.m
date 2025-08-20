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

        %% set up feasibility restoration phase 
       
        base_x = sparsity(nlsos.x);

        I = true(length(nlsos.g),1);
        % for idx = 1:length(nlsos.g)
        %     I(idx) = ~is_linear(nlsos.g(idx),nlsos.x);
        % end

        % get gram half-basis for nonlinear constraints
        [~,~,z] = grambasis(nlsos.g,I);
        
        % build unit vectors
        base_s0 = gramunit(z);
        
        r  = casos.PS.sym('r',sum(I));
   
        s0 = casos.PD(base_s0);

        nlsos_feas.x = [r;nlsos.x];
        
        nlsos_feas.g = nlsos.g(I) + r.*s0;
        
        x_R   = casos.PS.sym('x_R',base_x);

        if strcmp(obj.opts.feasibility_restoration,'Simple')
            % we simply mininimize on the constraint manifold
            Phi = [];
            lambda = [];
        elseif strcmp(obj.opts.feasibility_restoration,'Regularize')
            % check how far we are from the original problem/solution
            e   = nlsos_feas.x(sum(I)+1:end) - nlsos.x;
            Phi = 1/2*dot(e,e);
            lambda   = casos.PS.sym('l');
        elseif   strcmp(obj.opts.feasibility_restoration,'Cost')
            Phi = 1/2*nlsos.f;
            lambda   = casos.PS.sym('l');
        end

        nlsos_feas.f = sum(r) + lambda*Phi ;
        
        nlsos_feas.p = [nlsos.p; x_R;lambda];
        
        obj.FeasRes_para.n_r = sum(I);
        obj.FeasRes_para.length_dualOut = length(nlsos_feas.x)-length(nlsos_feas.g);
       
        sosopt.sossol         = obj.opts.sossol;
        sosopt.sossol_options = obj.opts.sossol_options;
        sosopt.Kx.lin         = length(nlsos_feas.x);
        sosopt.Kc.sos         = length(nlsos_feas.g);
        sosopt.error_on_fail  = false;
        sosopt.verbose        = 1;
        sosopt.max_iter       = 100;
        obj.feas_res_solver   =  casos.package.solvers.FeasibilityRestoration('feasRes',nlsos_feas,sosopt);
         
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




