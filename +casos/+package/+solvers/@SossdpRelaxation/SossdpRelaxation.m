classdef SossdpRelaxation < casos.package.solvers.SosoptCommon & matlab.mixin.Copyable
% Solve sum-of-squares problems by relaxation to SDP.

properties (Access=private)
    sdpsolver;

    gram2sos;
end

properties (Constant,Access=protected)
    sossdp_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'sdpsol_options', 'Options to be passed to the SDP solver.'}
    ];

    allow_eval_on_basis = true;
end

properties (SetAccess=private)
    class_name = 'SossdpRelaxation';
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SossdpRelaxation.sossdp_options;
    end
end

methods (Access=protected)
    % internal evaluation
    out = eval_on_basis(obj,in);
end

methods
    function obj = SossdpRelaxation(name,solver,sos,varargin)
        obj@casos.package.solvers.SosoptCommon(name,sos,varargin{:});

        % default options
        if ~isfield(obj.opts,'sdpsol_options'), obj.opts.sdpsol_options = struct; end
        % pass options to sdpsol
        if ~isfield(obj.opts.sdpsol_options,'error_on_fail')
            obj.opts.sdpsol_options.error_on_fail = obj.opts.error_on_fail;
        end

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
            % feasibility problem
            sos.f = casos.PS(0);
        else
            sos.f = casos.PS(sos.f);
        end
        % constraints
        if ~isfield(sos,'g')
            % unconstrained optimization
            sos.g = casos.PS;
        else
            sos.g = casos.PS(sos.g);
        end

        % build SDP problem
        buildproblem(obj,solver,sos);
    end

    function s = get_stats(obj)
        % Return stats.
        s = obj.sdpsolver.stats;
    end
end

methods (Access={?casos.package.functions.FunctionCommon, ?casos.package.functions.FunctionWrapper})
    function S = substitute(obj,idx,expr_in,expr_out)
        % Substitute a variable for expr_in -> expr_out.
        if ischar(idx)
            % map variable name to index
            idx = get_index_in(obj,idx);
        end

        % only parameter subsitution supported
        assert(idx == 1,'Subsitution of input %d not allowed.',idx)

        % project to basis
        [Qin,Zin] = poly2basis(expr_in);
        Qout = poly2basis(expr_out,obj.sparsity_p);

        % substitute
        S = copy(obj);
        S.sdpsolver = substitute(obj.sdpsolver,idx,Qin,Qout);
        % store new basis
        S.sparsity_p = Zin;
    end
end

end
